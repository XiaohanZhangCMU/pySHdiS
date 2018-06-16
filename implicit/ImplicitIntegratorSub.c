/**************************************************************************
 *
 *      Module:      ImplicitIntegratorSub.c
 *      Description: Used with the SubcycleIntegratorForceB function (force based subcycling).
 *                   Implements a numerical integrator using
 *                   the implicitly-solved trapezoid integration method.
 *					 This couples the trapezoid method with the Newton-Raphson
 *					 method for nonlinear system solution. The difference
 *					 between this function and ImplicitIntegrator is that it 
 *					 only time integrates the set of forces defined by reqType. 
 *
 ***************************************************************************/
#ifdef _IMPLICIT

#include "Home.h"
#include "sys/stat.h"
#include "sys/types.h"
#include "Implicit.h"
#include "Comm.h"
#include "cs.h"

#ifdef _LIS
#include "lis.h"
#endif

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

#define DIR_TIMESTEP_ERROR "timestep_error"
#define DIR_NODAL_TIMESTEP "timestep"

void   SortNodes(Home_t *home, real8 *maxsep);
void UnpackQMatrices(Node_t *node, real8 Q[3][3], real8 Qinv[3][3]);
void ConstrainJacobian(real8 Qseg[3][3], real8 Qinvseg[3][3], real8 Q1[3][3],
						real8 Qinv1[3][3], real8 Q2[3][3], real8 Qinv2[3][3],
						real8 J0[3][3], int ndof1, int ndof2, int locs[4],
						real8 J01a[3][3], real8 J01b[3][3], real8 J02a[3][3],
						real8 J02b[3][3]);
void AddJterms(Param_t *param, Node_t *node1, Node_t *node2, real8 subJ[3][3], cs *Jmob);



/*------------------------------------------------------------------------
 *
 *      Function:    PreserveNodalData
 *      Description: Both old and new values for certain nodal
 *                   data items are required during timestep
 *                   integration and calculating plastic strain.
 *                   This function copies the values for specified
 *                   items into the appropriate variables.
 *                   
 *-----------------------------------------------------------------------*/
static void PreserveNodalData(Home_t *home)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
 
/*
 *          Save previous nodal position.
 */
			node->oldx = node->x;
			node->oldy = node->y;
			node->oldz = node->z;

/*
 *          Make a copy of the current velocity.
 */
			node->currvX = node->vX;    node->oldvX = node->vX;
			node->currvY = node->vY;    node->oldvY = node->vY;
			node->currvZ = node->vZ;    node->oldvZ = node->vZ;
        }

        node = home->ghostNodeQ;

        while (node != (Node_t *)NULL) {

			node->oldx = node->x;
			node->oldy = node->y;
			node->oldz = node->z;

			node->currvX = node->vX;    node->oldvX = node->vX;
			node->currvY = node->vY;    node->oldvY = node->vY;
			node->currvZ = node->vZ;    node->oldvZ = node->vZ;

            node = node->next;
        }

        return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    UpdateJacobianSub
 *      Description: Update the mobility law Jacobian, Jmob. Jmob is stored
 *		using the coordinate (COO) format whereby every element is defined by
 *		three vectors - one defining the row, Jrow, one defining the comlumn
 *		Jcol, and one defining the value, Jval. For details see the lis
 *		user's manual. 
 *
 *
 *-----------------------------------------------------------------------*/
void UpdateJacobianSub(Home_t *home, cs *Jmob , 
                       Segm_t *SegList , int SegListSize)
{		
		int 	i, j, k, l, i1, i2, ndof1, ndof2, locs[4];
		real8	Lsum, Lsum1, Lsum2, bmag2, delta;
		real8	xtmp, ytmp, ztmp;
		real8	xdiff, ydiff, zdiff;
		real8	xdirseg[3], ydirseg[3], zdirseg[3];	
		real8	Qseg[3][3], Qinvseg[3][3];
		real8	Q1[3][3], Qinv1[3][3];	
		real8	Q2[3][3], Qinv2[3][3];
		real8	vectmp[3];
		real8	stressMat[3][3], stressMatloc[3][3];
		real8	tmpMat[3][3], tmpMat2[3][3], J0[3][3];
		real8	J01a[3][3], J01b[3][3], J02a[3][3], J02b[3][3];
		real8	subJ[3][3], JL[6][6], Japp[6][6], Jcore[6][6], Jseg[6][6];
		real8	r1_loc[3], r2_loc[3], dr[3], b_loc[3], dr_norm[3]; 
		real8	f1_loc[3], f2_loc[3], fapp[3];
		real8	L, coreFactor, s, B, Ecore;

		Node_t	*node, *nbr, *nbrtmp, *node1, *node2;
		Param_t *param;

		real8 	tol = 1e-10;
		int  	N = home->newNodeKeyPtr;
		param = home->param;

		real8 	*Lsumvec = malloc((N+1) * sizeof(real8));

/*
 *		We are updating the whole Jacobian, so flag every node.
 * 		Also calculate the segment length sum for all nodes that will be updated.
 */
		for (i=0; i < N; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node->flags |= UPDATE_NODE_J;
			Lsumvec[node->myTag.index] = 0.0;
		}
		
		for (i = 0; i < SegListSize; i++) {
			if (SegList[i].flag == 0) continue;
			
			node   = SegList[i].seg->node1;
            nbrtmp = SegList[i].seg->node2;
			
			xtmp = nbrtmp->x;
			ytmp = nbrtmp->y;
			ztmp = nbrtmp->z;
			PBCPOSITION(param, node->x, node->y, node->z, &xtmp, &ytmp, &ztmp);
			xdiff = node->x - xtmp;
			ydiff = node->y - ytmp;
			zdiff = node->z - ztmp;
			Lsum = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
			
			Lsumvec[node->myTag.index] += Lsum;
			Lsumvec[nbrtmp->myTag.index] += Lsum;
		}

/*
 *		Loop through the nodes and update the Jacobian. Nodes that have
 *		undergone topological changes (remeshing or collisions/splitting)
 *		need to be updated. When changing a Jmob element, we must first check
 * 		to see if it is already defined before appending the Jmob arrays to 
 *		add its value. 
 */

/* 
 *		We will need the current stress tensor later, so unpack it here.
 */

		stressMat[0][0] = param->appliedStress[0];
		stressMat[1][1] = param->appliedStress[1];
		stressMat[2][2] = param->appliedStress[2];
		stressMat[0][1] = param->appliedStress[5];
		stressMat[0][2] = param->appliedStress[4];
		stressMat[1][0] = param->appliedStress[5];
		stressMat[1][2] = param->appliedStress[3];
		stressMat[2][0] = param->appliedStress[4];
		stressMat[2][1] = param->appliedStress[3];

/*
 *		If a node and its neighbor are flagged, we add the contributions for 
 *		that segment. If only the main node is flagged, we add only its
 *		contribution.
 */
		for (i = 0; i < SegListSize; i++) {
			if (SegList[i].flag == 0) continue;
			
			node = SegList[i].seg->node1;
            nbr  = SegList[i].seg->node2;
			j = GetArmID(home, node, nbr);


			int tag1 = node->myTag.index;
			int tag2 = nbr->myTag.index;

/*
*			Calculate the lengths of all segments attached to both nodes.
*/
			Lsum1 = Lsumvec[node->myTag.index];
			Lsum2 = Lsumvec[nbr->myTag.index];

/*
 *			Calculate the local coordinate systems for both nodes
 *			and the segment.
 *			If a node only has one glide plane then the coordinate
 *			directions are as follows:
 *				x = Burgers vector direction
 *				z = glide plane normal
 *				y = cross product of z with x
 *			This is also the convention for the segment.
 *			If a node has more than one glide plane, then define
 *			the coordinate directions as follows:
 *				x = cross product of any two glide plane normals
 *				y = normal to any glide plane
 *				z = cross product of x with y
 *			Also calculate the resulting transformation matrices
 *			and store them for later use.
 */
			xdirseg[0] = node->burgX[j];
			xdirseg[1] = node->burgY[j];
			xdirseg[2] = node->burgZ[j];
			NormalizeVec(xdirseg);
			zdirseg[0] = node->nx[j];
			zdirseg[1] = node->ny[j];
			zdirseg[2] = node->nz[j];
			NormalizedCrossVector(zdirseg, xdirseg, ydirseg);
			Qseg[0][0] = xdirseg[0];
			Qseg[1][0] = xdirseg[1];
			Qseg[2][0] = xdirseg[2];
			Qseg[0][1] = ydirseg[0];
			Qseg[1][1] = ydirseg[1];
			Qseg[2][1] = ydirseg[2];
			Qseg[0][2] = zdirseg[0];
			Qseg[1][2] = zdirseg[1];
			Qseg[2][2] = zdirseg[2];

			Matrix33Transpose(Qseg, Qinvseg);

/*
 *			Transform everything into local coordinates for Jacobiancalculation.
 */
			UnpackQMatrices(node, Q1, Qinv1);
			UnpackQMatrices(nbr, Q2, Qinv2);

			vectmp[0] = node->x;
			vectmp[1] = node->y;
			vectmp[2] = node->z;
			Matrix33Vector3Multiply(Qinvseg, vectmp, r1_loc);
			xtmp = nbr->x;
			ytmp = nbr->y;
			ztmp = nbr->z;
			PBCPOSITION(param, node->x, node->y, node->z, &xtmp, &ytmp, &ztmp);
			vectmp[0] = xtmp;
			vectmp[1] = ytmp;
			vectmp[2] = ztmp;
			Matrix33Vector3Multiply(Qinvseg, vectmp, r2_loc);

			dr[0] = r2_loc[0] - r1_loc[0];
			dr[1] = r2_loc[1] - r1_loc[1];
			dr[2] = r2_loc[2] - r1_loc[2];
			L = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
			dr_norm[0] = dr[0]/L;
			dr_norm[1] = dr[1]/L;
			dr_norm[2] = dr[2]/L;

			vectmp[0] = node->burgX[j];
			vectmp[1] = node->burgY[j];
			vectmp[2] = node->burgZ[j];
			Matrix33Vector3Multiply(Qinvseg, vectmp, b_loc);

			Matrix33Mult33(Qinvseg, stressMat, tmpMat);
			Matrix33Mult33(tmpMat, Qseg, stressMatloc);
			
			vectmp[0] = node->fX;
			vectmp[1] = node->fY;
			vectmp[2] = node->fZ;
			Matrix33Vector3Multiply(Qinvseg, vectmp, f1_loc);

			vectmp[0] = nbr->fX;
			vectmp[1] = nbr->fY;
			vectmp[2] = nbr->fZ;
			Matrix33Vector3Multiply(Qinvseg, vectmp, f2_loc);

/*
 *			Calculate each of the partial Jacobians: core energy
 *			(line tension), length change, and applied stress. For
 *			each partial Jacobian, there are four sub-Jacobians 
 *			corresponding to 1-1, 1-2, 2-1, and 2-2 interactions.
 *			If either node is constrained to move in its local 
 *			x direction only (xonly# flag), then illegal terms need
 *			to be zeroed out using the ConstrainJacobian function.
 */
			B = 1/param->MobEdge;
			ndof1 = node->ndof;
			ndof2 = nbr->ndof;

/*
 *			Core energy partial Jacobian
 */
			if (param->elasticinteraction) {
				Ecore = param->Ecore;
			} else {
				Ecore = 0.5*param->TensionFactor*param->shearModulus;
			}

			bmag2 =	(node->burgX[j]*node->burgX[j] + node->burgY[j]*node->burgY[j] + 
							node->burgZ[j]*node->burgZ[j]);
			coreFactor = Ecore*bmag2/pow(L,5)*
					((1+param->pois)*dr[0]*dr[0]+(1-2*param->pois)*dr[1]*dr[1])/
					(1-param->pois);

			J0[0][0] = -coreFactor*dr[1]*dr[1];
			J0[0][1] =  coreFactor*dr[0]*dr[1];
			J0[1][0] =  coreFactor*dr[0]*dr[1];
			J0[1][1] = -coreFactor*dr[0]*dr[0];
			J0[0][2] = 0.0;
			J0[1][2] = 0.0;
			J0[2][0] = 0.0;
			J0[2][1] = 0.0;
			J0[2][2] = 0.0;

			locs[0] = 1; locs[1] = 1; locs[2] = 1; locs[3] = 1;

			ConstrainJacobian(Qseg, Qinvseg, Q1, Qinv1, Q2, Qinv2,
				J0, ndof1, ndof2, locs, J01a, J01b, J02a, J02b);
			
			for (i1=0; i1<3; i1++) {
				for (i2=0; i2<3; i2++) {
					Jcore[i1][i2]     =  J01a[i1][i2] / (B*Lsum1/2);
					Jcore[i1][i2+3]   = -J01b[i1][i2] / (B*Lsum1/2);
					Jcore[i1+3][i2]   = -J02a[i1][i2] / (B*Lsum2/2);
					Jcore[i1+3][i2+3] =  J02b[i1][i2] / (B*Lsum2/2);
				}
			} 

/*
 *			Applied stress partial Jacobian
 */
			Matrix33Vector3Multiply(stressMatloc, b_loc, vectmp);
			cross(vectmp, dr_norm, fapp);
			fapp[2] = 0;
			cross(dr_norm, fapp, vectmp);
			J0[0][0] = 0; J0[0][2] = 0; 
			J0[1][1] = 0; J0[1][2] = 0;
			J0[2][0] = 0; J0[2][1] = 0; J0[2][2] = 0;
			J0[0][1] = -0.5*vectmp[2];
			J0[1][0] =  0.5*vectmp[2];

			locs[0] = 0; locs[1] = 1; locs[2] = 0; locs[3] = 1;

			ConstrainJacobian(Qseg, Qinvseg, Q1, Qinv1, Q2, Qinv2,
				J0, ndof1, ndof2, locs, J01a, J01b, J02a, J02b);
			
			for (i1=0; i1<3; i1++) {
				for (i2=0; i2<3; i2++) {
					Japp[i1][i2]     =  J01a[i1][i2] / (B*Lsum1/2);
					Japp[i1][i2+3]   =  J01b[i1][i2] / (B*Lsum1/2);
					Japp[i1+3][i2]   = -J02a[i1][i2] / (B*Lsum2/2);
					Japp[i1+3][i2+3] =  J02b[i1][i2] / (B*Lsum2/2);
				}
			}

/*
 * 			Length change partial Jacobian
 */
			Matrix31Vector3Mult(f1_loc, dr_norm, tmpMat); 
			for (i1=0; i1<3; i1++) {
				for (i2=0; i2<3; i2++) {
					J0[i1][i2] = -tmpMat[i1][i2]/Lsum1;
				}
			}

			locs[0] = 1; locs[1] = 1; locs[2] = 0; locs[3] = 0;

			ConstrainJacobian(Qseg, Qinvseg, Q1, Qinv1, Q2, Qinv2,
				J0, ndof1, ndof2, locs, J01a, J01b, tmpMat, tmpMat2);
							
			
			Matrix31Vector3Mult(f2_loc, dr_norm, tmpMat); 
			for (i1=0; i1<3; i1++) {
				for (i2=0; i2<3; i2++) {
					J0[i1][i2] = tmpMat[i1][i2]/Lsum2;
				}
			}

			locs[0] = 0; locs[1] = 0; locs[2] = 1; locs[3] = 1;

			ConstrainJacobian(Qseg, Qinvseg, Q1, Qinv1, Q2, Qinv2,
				J0, ndof1, ndof2, locs, tmpMat, tmpMat2, J02a, J02b);
			
			for (i1=0; i1<3; i1++) {
				for (i2=0; i2<3; i2++) {
					JL[i1][i2]     = -J01a[i1][i2] / (B*Lsum1/2);
					JL[i1][i2+3]   =  J01b[i1][i2] / (B*Lsum1/2);
					JL[i1+3][i2]   =  J02a[i1][i2] / (B*Lsum2/2);
					JL[i1+3][i2+3] = -J02b[i1][i2] / (B*Lsum2/2);
					}
			}
/*
 *			Add them all together to get the segment Jacobian.
 */
			for (i1=0; i1<6; i1++) {
				for (i2=0; i2<6; i2++) {
					Jseg[i1][i2] = Jcore[i1][i2] + Japp[i1][i2] + JL[i1][i2];
				}
			}

/*
 *			Add these contributions to the global Jacobian. If node2 isn't flagged,
 *			then we don't add its terms.
 */
			subJ[0][0] = Jseg[0][0]; subJ[0][1] = Jseg[0][1]; subJ[0][2] = Jseg[0][2];
			subJ[1][0] = Jseg[1][0]; subJ[1][1] = Jseg[1][1]; subJ[1][2] = Jseg[1][2];
			subJ[2][0] = Jseg[2][0]; subJ[2][1] = Jseg[2][1]; subJ[2][2] = Jseg[2][2];
			AddJterms(param, node, node, subJ, Jmob);

			subJ[0][0] = Jseg[0][3]; subJ[0][1] = Jseg[0][4]; subJ[0][2] = Jseg[0][5];
			subJ[1][0] = Jseg[1][3]; subJ[1][1] = Jseg[1][4]; subJ[1][2] = Jseg[1][5];
			subJ[2][0] = Jseg[2][3]; subJ[2][1] = Jseg[2][4]; subJ[2][2] = Jseg[2][5];
			AddJterms(param, node, nbr, subJ, Jmob);

			subJ[0][0] = Jseg[3][0]; subJ[0][1] = Jseg[3][1]; subJ[0][2] = Jseg[3][2];
			subJ[1][0] = Jseg[4][0]; subJ[1][1] = Jseg[4][1]; subJ[1][2] = Jseg[4][2];
			subJ[2][0] = Jseg[5][0]; subJ[2][1] = Jseg[5][1]; subJ[2][2] = Jseg[5][2];
			AddJterms(param, nbr, node, subJ, Jmob);

			if ((nbr->flags & UPDATE_NODE_J)!=0) {
				subJ[0][0] = Jseg[3][3]; subJ[0][1] = Jseg[3][4]; subJ[0][2] = Jseg[3][5];
				subJ[1][0] = Jseg[4][3]; subJ[1][1] = Jseg[4][4]; subJ[1][2] = Jseg[4][5];
				subJ[2][0] = Jseg[5][3]; subJ[2][1] = Jseg[5][4]; subJ[2][2] = Jseg[5][5];
				AddJterms(param, nbr, nbr, subJ, Jmob);
			}
		} //loop over all nodes

/*
 * 		Unflag all nodes.
 */
		for (i=0; i < N; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            node->flags &= ~UPDATE_NODE_J;
		}

/*
 * 		Clean up memory.
 */
		free(Lsumvec);

		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    ImplicitIntegratorSub
 *      Description: Implements a numerical timestep integrator using
 *                   the implicitly-solved Trapezoid integration method.
  *					 Only integrates group 0 set of forces.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
static int Implicit_Step(Home_t *home , int reqType , real8 newDT , int *incrDelta)
{
        int      i, j, convergent, maxIterations, tmp;
        int 	 ind1, ind2, ind3;
        int      iter, globalIterError, mobIterError;
        int      dumpErrorData = 1, doAll = 1;
		int 	 errInc, count, Ndof;
		real8    errMax, globalErrMax, oldglobalErrMax, drn, relerrMax;
        real8    oldx, oldy, oldz, globalRelErrMax;
        real8    localVals[3], globalVals[3];
		real8 	 Q[3][3], Qinv[3][3];
		real8	 rtmp[3], rout[3], vtmp[3], vout[3];
		real8	 thisx, thisy, thisz;
		real8 	 errdX, errdY, errdZ, f1, f2, f3;
		real8	 dx, dy, dz, xtmp, ytmp, ztmp;
        Node_t   *node;
        Param_t  *param;
		Subcyc_t *subcyc;

        param  = home->param;
		subcyc = home->subcyc;
		int N  = home->newNodeKeyPtr;

#ifndef _LIS
		Fatal("Linear Algebra Library (-D_LIS) is required for implicit integrator!");
                return (0);
#else
		LIS_VECTOR 	negf_lis, dr_lis;
		LIS_MATRIX 	J_lis, Jmob_lis;
		LIS_SOLVER 	solver;

/*
 * 		Allocate memory for arrays.
 */
		cs      *J, *Jmobtrp, *Jmob, *Eye, *tmp_cs;
		real8	*f, *r0, *v0, *Jval;
		int 	*Jrow, *Jcol;
		tmp = 3*(N+1);
		f   = (real8 *)malloc(tmp * sizeof(real8));
		r0  = (real8 *)malloc(tmp * sizeof(real8));
		v0  = (real8 *)malloc(tmp * sizeof(real8));

/*
 *		Set the initial guess to the current nodal positions and velocities
 */
		for (i=0; i < N; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if (node->ndof==0) continue;

			node->x = node->oldx;   node->vX = node->currvX;
			node->y = node->oldy;   node->vY = node->currvY;
			node->z = node->oldz;   node->vZ = node->currvZ;
		}
		
		mobIterError  = 0;
		errInc        = 0;
		globalVals[0] = 0.0;
		globalVals[1] = 0.0;
		globalVals[2] = 0.0;

/*
 *		Calculate nodal forces and velocities
 */
    /*	NodeForce(home, reqType);
        mobIterError = CalcNodeVelocities(home, 0, doAll);
        CommSendVelocity(home);
		if (mobIterError) return(0);*/

/*
 *      Loop until we converge on a time step.  
 */
        convergent    = 0;
        maxIterations = 25;

/*
 * 		Calculate coordinate systems and transformation matrices for every
 * 		node.
 */		
		for (i=0; i < N; i++) {
	        node = home->nodeKeys[i];
	        if (node == (Node_t *)NULL) continue;
			node->flags |= UPDATE_NODE_J;
		}
		Ndof = 0;
		AssignCoordSys(home, &Ndof, -1);
		
/*
 *		Calculate the Jacobian
 */ 
		Jmobtrp = cs_spalloc(Ndof, Ndof, 2*Ndof, 1, 1);
		if (reqType == GROUP0) {
			UpdateJacobianSub(home,Jmobtrp, subcyc->SegListG0, subcyc->SegListG0_cnt);
		} else if (reqType == GROUP1) {
			UpdateJacobianSub(home,Jmobtrp, subcyc->SegListG1, subcyc->SegListG1_cnt);
		} else if (reqType >= GROUP2) {
			UpdateJacobianSub(home,Jmobtrp, NULL     , 0            );
		}
		int 	changedJ;
	//	UpdateJacobian(home,1,&changedJ,Jmobtrp); // Jmobtrp gets updated here. 1 means update all nodes. changedJ is useless.
		
		Jmob = cs_compress(Jmobtrp);
		cs_spfree(Jmobtrp);
		cs_dupl(Jmob);
		
		tmp_cs = cs_spalloc(Ndof, Ndof, Ndof, 1, 1) ;    /* create triplet identity matrix */
		for (i = 0 ; i < Ndof ; i++) cs_entry(tmp_cs, i, i, 1) ;
		Eye = cs_compress(tmp_cs);
		cs_spfree(tmp_cs);

/*
 *		Calculate the iteration matrix. Diagonal entries first, then off-
 *		diagonal entries.
 */		
		lis_matrix_create(0,&J_lis);
		lis_matrix_set_size(J_lis,0,Ndof);

		J = cs_add(Eye, Jmob, 1, -newDT/2);

		int column = 0;
		for (j=0; j < J->nzmax; j++) {
			if (j == J->p[column+1]) column++;
			
			lis_matrix_set_value(0,J->i[j], column, J->x[j], J_lis);
		}
		lis_matrix_assemble(J_lis);

/*
 *		Set the initial guess to the current nodal positions.
 *		Convert the old position and velocity vectors to local coordinates
 */
		for (i=0; i < N; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if ( node->ndof==0) continue;

			UnpackQMatrices(node, Q, Qinv);

			rtmp[0] = node->oldx;   vtmp[0] = node->currvX;
			rtmp[1] = node->oldy;   vtmp[1] = node->currvY;
			rtmp[2] = node->oldz;   vtmp[2] = node->currvZ;					

			Matrix33Vector3Multiply(Qinv, rtmp, rout);
			Matrix33Vector3Multiply(Qinv, vtmp, vout);
			
			r0[i    ] = rout[0];   v0[i    ] = vout[0];
			r0[i+N  ] = rout[1];   v0[i+N  ] = vout[1];
			r0[i+2*N] = rout[2];   v0[i+2*N] = vout[2];
		}
		
		for (iter = 0; iter < maxIterations; iter++) {
/*
 *			Recalculate nodal force and velocity at the new positions
 */
			if (iter>0) { //No need to update velocities on 1st iteration
				NodeForce(home, reqType);
				mobIterError = CalcNodeVelocities(home, 0, doAll);
				CommSendVelocity(home);
				if (mobIterError) return(0);
			}

/*
 *			Calculate the residual vector and maximum error in local
 *			coordinates.
 */
			errMax    = 0.0;
			relerrMax = 0.0;
			
			for ( i=0 ; i < N ; i++) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				if ( node->ndof == 0) node->subgroup = 0;
				if ( node->ndof == 0) continue;

				UnpackQMatrices(node, Q, Qinv);

				thisx = node->x;
				thisy = node->y;
				thisz = node->z;
				PBCPOSITION(param, node->oldx, node->oldy, node->oldz,
				                   &thisx    , &thisy    , &thisz   );

				rtmp[0] = thisx;
				rtmp[1] = thisy;
				rtmp[2] = thisz;					
				vtmp[0] = node->vX;
				vtmp[1] = node->vY;
				vtmp[2] = node->vZ;

				Matrix33Vector3Multiply(Qinv, rtmp, rout);
				Matrix33Vector3Multiply(Qinv, vtmp, vout);

				if (node->dofx != 0) {
					f1 = rout[0] - r0[i] - newDT/2*(vout[0] + v0[i]);
					f[node->dofx-1] = f1;
				} else {
					f1 = 0;
				}
				if (node->dofy != 0) {
					f2 = rout[1] - r0[i+N] - newDT/2*(vout[1] + v0[i+N]);
					f[node->dofy-1] = f2;
				} else {
					f2 = 0;
				}
				if (node->dofz != 0) {
					f3 = rout[2] - r0[i+2*N] - newDT/2*(vout[2] + v0[i+2*N]);
					f[node->dofz-1] = f3;
				} else {
					f3 = 0;
				}
				
				real8 thiserr = sqrt( f1*f1 + f2*f2 + f3*f3);
				errMax = MAX(errMax, thiserr);
				
				//Calculate relative error
				errdX = thisx - node->oldx;
				errdY = thisy - node->oldy;
				errdZ = thisz - node->oldz;
				drn   = sqrt( errdX*errdX + errdY*errdY + errdZ*errdZ );
				
				if (thiserr > param->rTolth) {
					if (drn > param->rTolth / param->rTolrel) {
						relerrMax = MAX( relerrMax , thiserr/drn);
					} else {
						relerrMax = 2.0 * param->rTolrel;
					}
				}
				
				if ((thiserr < param->rTol) && (thiserr < param->rTolth || thiserr/drn < param->rTolrel) && iter>0) {
					node->subgroup = 0;
				} else {
					node->subgroup = 1;
				}
			}

/*
 *			Need to find largest errMax from among all domains.
 */
#if PARALLEL
			localVals[0] = errMax;
			localVals[1] = (real8)mobIterError;
			localVals[2] = relerrMax;

			MPI_Allreduce(localVals, globalVals, 3, MPI_DOUBLE, MPI_MAX, 
			              MPI_COMM_WORLD);

			globalErrMax    = globalVals[0];
			globalIterError = globalVals[1];
			globalRelErrMax = globalVals[2];
#else
			globalErrMax    = errMax;
			globalIterError = mobIterError;
			globalRelErrMax = relerrMax;
#endif

/*
 *			If any domain encountered an error iterating inside
 *			the mobility function, just go right to cutting the
 *			timestep and trying again.
 */
			if (globalIterError) return(0);

/*
 * 			If the error is within the tolerance, we've reached
 *			convergence so we can accept this deltaT.  Otherwise
 *			continue iterating, unless max number of iterations has reached.
 *			If it has, exit the for loop and cut the time step.  
 *			Note: we need to reposition both local nodes and ghost nodes!
 */
			if (globalErrMax < param->rTol && globalRelErrMax < param->rTolrel && iter > 0) {
				if (iter>5) *incrDelta = 0;
				convergent = 1;
				break;
			} else {
				if (iter == maxIterations-1) continue;
				
				if (iter>0) {
					if (globalErrMax > oldglobalErrMax) {
						errInc++;
						if (errInc>=3) break;
					}
				}
				oldglobalErrMax = globalErrMax;

/*				Solve for the new nodal positions and velocities using the 
 *				Newton-Raphson method.
 */
				//Instantiate lis vectors and solver.
				lis_vector_create(0, &negf_lis);
				lis_vector_create(0, &dr_lis);

				lis_solver_create(&solver);
				lis_solver_set_option("-i bicg -p none",solver);
				lis_solver_set_option("-tol 1.0e-12",solver);

				lis_vector_set_size(negf_lis, Ndof, 0);
				lis_vector_set_size(dr_lis, Ndof, 0);

				//Set values of negf_lis vector.
				for (i=0; i < N; i++) {
		            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
					if ( node->ndof==0) continue;

					ind1 = node->dofx;
					ind2 = node->dofy;
					ind3 = node->dofz;

					if (ind1 != 0) lis_vector_set_value(0, ind1-1, -f[ind1-1], negf_lis);
					if (ind2 != 0) lis_vector_set_value(0, ind2-1, -f[ind2-1], negf_lis);
					if (ind3 != 0) lis_vector_set_value(0, ind3-1, -f[ind3-1], negf_lis);
				}

				//Solve the system
				lis_solve(J_lis,negf_lis,dr_lis,solver);

				//Extract the new nodal positions and rotate to global coords
				for (i=0; i < N; i++) {
					if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
					if ( node->ndof==0) continue;

					ind1 = node->dofx;
					ind2 = node->dofy;
					ind3 = node->dofz;

					if (ind1 != 0) lis_vector_get_value(dr_lis, ind1-1, &dx);
					else           dx = 0;
					
					if (ind2 != 0) lis_vector_get_value(dr_lis, ind2-1, &dy);
					else           dy = 0;
					
					if (ind3 != 0) lis_vector_get_value(dr_lis, ind3-1, &dz);
					else           dz = 0;

					UnpackQMatrices(node, Q, Qinv);

					rtmp[0] = dx;
					rtmp[1] = dy;
					rtmp[2] = dz;					

					Matrix33Vector3Multiply(Q, rtmp, rout);

					xtmp = node->x + rout[0];
					ytmp = node->y + rout[1];
					ztmp = node->z + rout[2];

					FoldBox(param,&xtmp,&ytmp,&ztmp);

					node->x = xtmp;
					node->y = ytmp;
					node->z = ztmp;					
				}

				//Kill lis vectors.
				lis_solver_destroy(solver);
				lis_vector_destroy(negf_lis);
				lis_vector_destroy(dr_lis);

			}  /* not convergent */
		}  /* for (iter = 0; ...) */
            
/*
 *		Need to kill the lis object because we will reinstantiate it.
 */
		lis_matrix_destroy(J_lis);

/*
 *		Clean up memory.
 */
		free(r0);
		free(v0);
		free(f);
		cs_spfree(Jmob);
		cs_spfree(J) ;
		cs_spfree(Eye); 
        
		return(convergent);
#endif
}


/*------------------------------------------------------------------------
 *
 *      Function:     cutTimeStepSize
 *      Description:  This function cuts the delta T by a configured 
 *                    factor, if we have not reached convergence
 *
 *-----------------------------------------------------------------------*/
static void cutTimeStepSize(Home_t *home , real8 *newDT , int *incrDelta)
{
		Param_t *param;
		param = home->param;

		*newDT         *=  param->dtDecrementFact;
		*incrDelta      =  0;
		param->deltaTT  = *newDT;

		if ((*newDT < 1.0e-20) && (home->myDomain == 0)) {
			Fatal("TrapezoidIntegrator(): Timestep has dropped below\n"
				  "minimal threshold to %e.  Aborting!", *newDT);
		}
}


/*------------------------------------------------------------------------
 *
 *      Function:    moveFromG0toG4
 *      Description: Moves the interactions of the glagged nodes
 *                   from group0 to group4
 *
 *-----------------------------------------------------------------------*/
void moveFromG0toG4(Home_t *home , int *incrDelta)
{
		int       i , mobIterError , doAll = 1;
		real8     rg9s;
		Node_t   *node1 , *node2 , *node3 , *node4 , *node;
		Param_t  *param ;
		Subcyc_t *subcyc;
		
	//	*incrDelta    =  0;
		param  = home->param;
		subcyc = home->subcyc;
		rg9s   = param->rg4 * param->rg4;
		
		for (i = 0; i < subcyc->SegListG0_cnt ; i++) {
			if (subcyc->SegListG0[i].flag == 0) continue;
			
			node1 = subcyc->SegListG0[i].seg->node1;
			node2 = subcyc->SegListG0[i].seg->node2;
			
			if ((node1->subgroup == 1) || (node2->subgroup == 1)) {
				subcyc->SegListG0[i].flag = 0;
				
				subcyc->SegListG1[subcyc->SegListG1_cnt].seg  = subcyc->SegListG0[i].seg;
				subcyc->SegListG1[subcyc->SegListG1_cnt].flag = 1;
				subcyc->SegListG1_cnt++;
				
			/*	if (node1->subgroup == 1) node1->ndof = 0;
				if (node2->subgroup == 1) node2->ndof = 0;*/
			}
		}
		
		for (i = 0; i < subcyc->SegSegListG0_cnt ; i++) {
			if (subcyc->SegSegListG0[i].flag == 0) continue;
			
			node1 = subcyc->SegSegListG0[i].seg1->node1;
			node2 = subcyc->SegSegListG0[i].seg1->node2;
			node3 = subcyc->SegSegListG0[i].seg2->node1;
			node4 = subcyc->SegSegListG0[i].seg2->node2;
			
			if ((node1->subgroup == 1) || (node2->subgroup == 1) ||
				(node3->subgroup == 1) || (node4->subgroup == 1)) {
				
				if (subcyc->SegSegListG0[i].dist2 > rg9s) continue;
				subcyc->SegSegListG0[i].flag = 0;
			
				if (subcyc->SegSegListG4 == NULL || subcyc->SegSegListG4_cnt >= subcyc->SegSegListG4_siz) {
					subcyc->SegSegListG4_siz += 100;
					subcyc->SegSegListG4 = (SegSeg_t *)realloc(subcyc->SegSegListG4, sizeof(SegSeg_t) *
													   subcyc->SegSegListG4_siz);
				}
				
				subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].seg1  = subcyc->SegSegListG0[i].seg1;
				subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].seg2  = subcyc->SegSegListG0[i].seg2;
				subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].flag  = 1;
				subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].setSeg1Forces = 1;
				subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].setSeg2Forces = 1;
				subcyc->SegSegListG4_cnt++;
				
			/*	if (node1->subgroup == 1) node1->ndof = 0;
				if (node2->subgroup == 1) node2->ndof = 0;
				if (node3->subgroup == 1) node3->ndof = 0;
				if (node4->subgroup == 1) node4->ndof = 0;*/
			}
		}
		
/*
 *		Since GROUP0 has changed we need to calculate the current forces 
 *		and velocities again
 */
		for ( i = 0 ; i < home->newNodeKeyPtr ; i++ ) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node->x = node->oldx;
			node->y = node->oldy;
			node->z = node->oldz;
		}
			
		NodeForce(home, GROUP0);
		mobIterError = CalcNodeVelocities(home, 0, doAll);
		CommSendVelocity(home);
		
		for ( i = 0 ; i < home->newNodeKeyPtr ; i++ ) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node->currvX = node->vX;    node->oldvX = node->vX;
			node->currvY = node->vY;    node->oldvY = node->vY;
			node->currvZ = node->vZ;    node->oldvZ = node->vZ;
		}
		
		return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    ImplicitIntegratorSub
 *      Description: Implements a numerical timestep integrator using
 *                   the implicitly-solved Trapezoid integration method.
 *					 Integrates the group defined by reqType.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
void ImplicitIntegratorSub(Home_t *home , int reqType)
{
		int     i, convergent, incrDelta, ioGroup = home->ioGroupNum;
		int     globalIterError, mobIterError;
		int     doAll = 1, N = home->newNodeKeyPtr;
		char    debug_data [50];
		real8   errMax, localVals[2], globalVals[2];
		real8   newDT , oldx  , oldy  , oldz     , drn;
		real8	errx  , erry  , errz  , thiserr  , factor;
		real8   errdX , errdY , errdZ , relerrMax, globalErrMax;
		real8   tmp1  , tmp2  , tmp3  , tmp4     , globalRelErrMax;
		Node_t  *node;
		Param_t *param;

		param = home->param;
		
/*
 *		Calculate the current velocities for the given group defined by reqType
 */
		NodeForce(home, reqType);
		mobIterError = CalcNodeVelocities(home, 0, doAll);
		CommSendVelocity(home);
		
/*
 *		Allocating the required arrays
 */
		int    *col = malloc(N * sizeof(int     ));
		real8  **X0 = malloc(N * sizeof(real8 * ));
		real8  **X1 = malloc(N * sizeof(real8 * ));
		
		for ( i = 0 ; i < N ; i++ ) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			X0[i] = malloc( 3 * sizeof(real8 ));
			X1[i] = malloc( 3 * sizeof(real8 ));
			
			X0[i][0] = node->x;
			X0[i][1] = node->y;
			X0[i][2] = node->z;
			col[i]   = 0;
		}

/*
 * 		Set the initial time step.
 */
		if      (reqType == GROUP0) newDT = MIN(param->maxDT, param->nextDT    );
		else if (reqType == GROUP1) newDT = MIN(param->maxDT, param->nextDTsub );
		else if (reqType == GROUP2) newDT = MIN(param->maxDT, param->nextDTsub2);
		else if (reqType == GROUP3) newDT = MIN(param->maxDT, param->nextDTsub3);
		else if (reqType == GROUP4) newDT = MIN(param->maxDT, param->nextDTsub4);
		
		if (newDT <= 0.0) newDT = param->maxDT;
		
		if      (reqType == GROUP0) param->deltaTT     = newDT;
		else if (reqType == GROUP1) param->deltaTTsub  = newDT;
		else if (reqType == GROUP2) param->deltaTTsub2 = newDT;
		else if (reqType == GROUP3) param->deltaTTsub3 = newDT;
		else if (reqType == GROUP4) param->deltaTTsub4 = newDT;

/*
 *      Preserve certain items of nodal data for later use.
 */
        PreserveNodalData(home);

/*
 *      Loop until we converge on a time step.  First step is to
 *      use the current positions and velocities and reposition the
 *      nodes to where they would be after the suggested delta time.
 */
        convergent = 0;
        incrDelta  = 1;
		
        int iTry1=0;
		int iTry2=0;
		
		while (!convergent) {
			errMax        = 0.0;
			relerrMax     = 0.0;
			globalVals[0] = 0.0;
			globalVals[1] = 0.0;
			iTry2++;
		
/*
 *			Taking one step of size dt and storing the positions in X1
 */
			convergent = Implicit_Step(home , reqType , newDT , &incrDelta);
			if (!convergent) {
				iTry1++;
			
				if (iTry2 <= param->nTry && reqType==GROUP0) {
					moveFromG0toG4 (home , &incrDelta);
					
				//	sprintf (debug_data, "first%d", home->cycle);
				//	WriteAtomEye(home, debug_data, ioGroup, 1, 1, 1);
				} else {
					cutTimeStepSize(home , &newDT , &incrDelta);
				}
				continue;
			}
		//	continue;
			
			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				X1[i][0] = node->x;
				X1[i][1] = node->y;
				X1[i][2] = node->z;
			}
			
/*
 *			Taking two implicit steps of size dt/2
 *			First take one step of size dt/2 starting form the old positions
 */
			convergent = Implicit_Step(home , reqType , newDT*0.5 , &incrDelta);
			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				node->oldx = node->x;   node->currvX = node->vX;
				node->oldy = node->y;   node->currvY = node->vY;
				node->oldz = node->z;   node->currvZ = node->vZ;
			}
			
/*
 *			Now take the second half step and evaluate the truncation error
 */
			convergent = Implicit_Step(home , reqType , newDT*0.5 , &incrDelta);

/*
 *			Restore the old positions and velocities
 */ 			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				node->oldx = X0[i][0];   node->currvX = node->oldvX;
				node->oldy = X0[i][1];   node->currvY = node->oldvY;
				node->oldz = X0[i][2];   node->currvZ = node->oldvZ;
			}
			
			if (!convergent) {
				iTry1++;
				
				if (iTry2 <= param->nTry && reqType == GROUP0) {
					moveFromG0toG4 (home , &incrDelta);
					
				//	sprintf (debug_data, "second%d", home->cycle);
				//	WriteAtomEye(home, debug_data, ioGroup, 1, 1, 1);
				} else {
					cutTimeStepSize(home , &newDT , &incrDelta);
				}
				continue;
			}
			
/*
 *			Calculate the truncation error on this node
 */			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				errx = node->x - X1[i][0];
				erry = node->y - X1[i][1];
				errz = node->z - X1[i][2];
				ZImage(param, &errx, &erry, &errz);
				
				thiserr = sqrt( errx*errx + erry*erry + errz*errz );
				errMax  = MAX ( errMax , thiserr);

/*
 *				Calculate the relative error on this node
 */
				errdX = node->x - X0[i][0];
				errdY = node->y - X0[i][1];
				errdZ = node->z - X0[i][2];
				ZImage(param, &errdX, &errdY, &errdZ);
				drn   = sqrt( errdX*errdX + errdY*errdY + errdZ*errdZ );

				if (thiserr > param->rTolth ) {
					if (drn > param->rTolth / param->rTolrel ) {
						relerrMax = MAX( relerrMax , thiserr/drn );
					} else {
						relerrMax = 2.0*param->rTolrel;
					}
				}
				
				if (iTry2-iTry1 <= param->nTry-1 && reqType == GROUP0) {
					if ((thiserr < param->rTol) && (thiserr < param->rTolth || thiserr/drn < param->rTolrel)) {
						node->subgroup = 0;
					} else {
						node->subgroup = 1;
						col[i] = 1;
					}
				}
			}
			
/*
 *			Need to find largest errMax from among all domains.
 */
#if PARALLEL
			localVals[0] = errMax;
			localVals[1] = relerrMax;

			MPI_Allreduce(localVals, globalVals, 2, MPI_DOUBLE, MPI_MAX, 
			              MPI_COMM_WORLD);

			globalErrMax    = globalVals[0];
			globalRelErrMax = globalVals[1];
#else
			globalErrMax    = errMax;
			globalRelErrMax = relerrMax;
#endif

/*
 *			If the error is within the tolerance, we've reached
 *			convergence so we can accept this deltaT.  Otherwise
 *			reposition the nodes and try again.  Note: we need to
 *			reposition both local nodes and ghost nodes!
 */
			if (globalErrMax < param->rTol && globalRelErrMax < param->rTolrel) {
				convergent = 1;
				break;
			} else {
				convergent = 0;
				
				if (iTry2-iTry1 <= param->nTry-1 && reqType == GROUP0) 
					moveFromG0toG4 (home , &incrDelta);
				else
					cutTimeStepSize(home , &newDT , &incrDelta);
			}
		}
		
		printf("iTry1 = %d  ,  iTry2 = %d\n", iTry1 , iTry2);
		
		int  plotAE = 0;
		for (i=0 ; i<N ; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node->subgroup = col[i];
			if ( col[i] > 0 ) plotAE = 1;
		}
		
		if (plotAE == 1) {
		//	sprintf (debug_data, "third%d", home->cycle);
		//	WriteAtomEye(home, debug_data, ioGroup, 1, 1, 1);
		}

/*
 *      Automatically increment timestep if convergence was reached
 *      on the very first iteration of the above loop.
 *
 *      If variable timestep adjustments are enabled, calculate an
 *      adjustment factor based on the maximum allowed timestep increment
 *      and the maximum error found above.  If variable timestep
 *      adjustments are not enabled, adjust the timestep by the
 *      maximum permitted factor.
 */
		
		if (reqType == GROUP0) {
			param->deltaTT     = newDT;
			param->realdt      = newDT;
			param->timeStart   = param->timeNow;
		} else if (reqType == GROUP1) {
			param->deltaTTsub  = newDT;
			param->realdtsub   = newDT;
		} else if (reqType == GROUP2) {
			param->deltaTTsub2 = newDT;
			param->realdtsub2  = newDT;
		} else if (reqType == GROUP3) {
			param->deltaTTsub3 = newDT;
			param->realdtsub3  = newDT;
		} else if (reqType == GROUP4) {
			param->deltaTTsub4 = newDT;
			param->realdtsub4  = newDT;
		}
		
		if (incrDelta) {
            if (param->dtVariableAdjustment) {
                
                tmp1 = pow(param->dtIncrementFact, param->dtExponent);
                tmp2 = globalErrMax/param->rTol;
                tmp3 = 1.0 / param->dtExponent;
                tmp4 = pow(1.0/(1.0+(tmp1-1.0)*tmp2), tmp3);
                factor = param->dtIncrementFact * tmp4;
				newDT  = MIN(param->maxDT, newDT*factor);
					
            } else {
				newDT = MIN(param->maxDT, newDT*param->dtIncrementFact);
            }	
        }
			
		if      (reqType == GROUP0) param->nextDT     = newDT;
		else if (reqType == GROUP1) param->nextDTsub  = newDT;
		else if (reqType == GROUP2) param->nextDTsub2 = newDT;
		else if (reqType == GROUP3) param->nextDTsub3 = newDT;
		else if (reqType == GROUP4) param->nextDTsub4 = newDT;
		
/*
 *		Freeing memory
 */ 
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			free( X0[i] );
			free( X1[i] );
		}
		free( X0 );
		free( X1 );
		free( col);
		
		return;
		
}
#endif //_IMPLICIT
