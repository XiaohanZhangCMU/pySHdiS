/**************************************************************************
 *
 *  Function    : Mobility_FCC_0
 *  Author      : Wei Cai, Seok-Woo Lee (updated 07/14/09)
 *  Description : Generic Mobility Law of FCC metals
 *                Each line has a glide plane from input
 *                and it never changes
 *                If the plane normal is not of {111} type, dislocation
 *                motion is constrained along line direction
 *                If node flag == 7, node velocity is zero
 *
 *  Returns:  0 on success
 *            1 if velcoity could not be determined
 *
 ***************************************************************************/

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif

#define _ENABLE_LINE_CONSTRAINT 1
#define _PRECOMPUTE_BMATRIX 0

/*-------------------------------------------------------------------------
 *
 *      Function:     Mobility_FCC_0
 *      Description:  Wrapper function to chose the appropriate function
 *                    to call to evaluate the nodal mobility
 * 
 *-----------------------------------------------------------------------*/
int Mobility_FCC_0(Home_t *home, Node_t *node)
{	
	int      mob;
	Param_t *param;
	
	param = home->param;
	
#ifdef _PRECOMPUTE_GLIDE
	if (param->preMobilityConstraints == 1 && node->mobMatrixUse == 1) {
		mob = Mobility_FCC_0_matrix(home, node);
	} else {
		mob = Mobility_FCC_0_original(home, node);
	}
#else		
	mob = Mobility_FCC_0_original(home, node);
#endif
	
	return mob;
}

#ifdef _GPU_SUBCYCLE
/*-------------------------------------------------------------------------
 *
 *      Function:     Mobility_FCC_0_matrix_GPU
 *      Description:  This function precomputes a nodal mobility matrix 
 *                    operator lumping the glide and line constraints
 *                    to evaluate the mobility from forces.
 * 
 *-----------------------------------------------------------------------*/
int Mobility_FCC_0_matrix_GPU(Home_t *home, Node_t *node, double mobMatrix[3][3])
{
		int numNonZeroLenSegs = 0;
		Param_t *param;
		real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
		int i, j, k, nc, nconstraint, nlc;
		real8 normX[100], normY[100], normZ[100];
		real8 normx[100], normy[100], normz[100];
		real8 lineX[100], lineY[100], lineZ[100];
		real8 a, b;
		real8 dx, dy, dz, lx, ly, lz, lr, LtimesB;
		real8 lcx, lcy, lcz, normdotlc, lc[3];
		Node_t *nbr;
		real8 MobScrew, MobEdge, Mob;
		real8 bx, by, bz, br, dangle;
		real8 nForce[3];

		param = home->param;

		MobScrew = param->MobScrew;
		MobEdge  = param->MobEdge;

		nc = node->numNbrs;
		
/*
 *  	Initialize mobility matrix operator
 */		
		//node->mobMatrixUse = 1;
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				mobMatrix[i][j] = 0.0;
  
/*
 *  	If node is 'pinned' in place by constraints, or the node has any arms 
 *  	with a burgers vector that has explicitly been set to be sessile (via
 *  	control file inputs), the node may not be moved so just zero the velocity
 *  	and return
 */
		if ((node->constraint == PINNED_NODE) || NodeHasSessileBurg(home, node)) {
			return 0;
		}

/*
 *  	It's possible this function was called for a node which had only zero-
 *  	length segments (during SplitSurfaceNodes() for example).  If that is
 *  	the case, just set the velocity to zero and return.
 */
		for (i = 0; i < nc; i++) {
			if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) continue;
			dx = node->x - nbr->x;
			dy = node->y - nbr->y;
			dz = node->z - nbr->z;
			if ((dx*dx + dy*dy + dz*dz) > 1.0e-12) {
				numNonZeroLenSegs++;
			}
		}

		if (numNonZeroLenSegs == 0) {
			return 0;
		}


/*
 *      Copy glide plane constraints and determine
 *      line constraints.
 */
		for(i=0;i<nc;i++) {
			
			normX[i] = normx[i] = node->nx[i];
			normY[i] = normy[i] = node->ny[i];
			normZ[i] = normz[i] = node->nz[i];

/*
 *      	If needed, rotate the glide plane normals from the
 *      	laboratory frame to the crystal frame.
 */
			if (param->useLabFrame) {
				real8 normTmp[3] = {normX[i], normY[i], normZ[i]};
				real8 normRot[3];

				Matrix33Vector3Multiply(home->rotMatrixInverse, normTmp, normRot);

				normX[i] = normRot[0]; normY[i] = normRot[1]; normZ[i] = normRot[2];
				normx[i] = normRot[0]; normy[i] = normRot[1]; normz[i] = normRot[2];
			}

			if ( (fabs(fabs(normX[i]) - fabs(normY[i])) > FFACTOR_NORMAL) ||
                 (fabs(fabs(normY[i]) - fabs(normZ[i])) > FFACTOR_NORMAL) ) {
				
				/* not {111} plane */
				if ((nbr=GetNeighborNode(home,node,i)) == (Node_t *)NULL) {
					Fatal("Neighbor not found at %s line %d\n",__FILE__,__LINE__);
				}
				
				lineX[i] = nbr->x - node->x;
				lineY[i] = nbr->y - node->y;
				lineZ[i] = nbr->z - node->z;
				ZImage (param, lineX+i, lineY+i, lineZ+i);

/*
 *          	If needed, rotate the line sense from the laboratory frame to
 *          	the crystal frame.
 */
				if (param->useLabFrame) {
					real8 lDir[3] = {lineX[i], lineY[i], lineZ[i]};
					real8 lDirRot[3];

					Matrix33Vector3Multiply(home->rotMatrixInverse, lDir, lDirRot);

					lineX[i] = lDirRot[0];
					lineY[i] = lDirRot[1];
					lineZ[i] = lDirRot[2];
				}
				
			} else { 
				/* no line constraint */
				lineX[i] = lineY[i] = lineZ[i] = 0;
			}
		}
    
/*
 * 		Normalize glide plane normal vectors 
 *		and lc line vectors
 */
		for(i=0;i<nc;i++) {
			
			a=sqrt(normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i]);
			b=sqrt(lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i]);

			if(a>0) {
				normX[i]/=a;
				normY[i]/=a;
				normZ[i]/=a;

				normx[i]/=a;
				normy[i]/=a;
				normz[i]/=a;
			}
			if(b>0) {
				lineX[i]/=b;
				lineY[i]/=b;
				lineZ[i]/=b;
			}
		}

/*
 * 		Find independent glide constraints
 */
		nconstraint = nc;
		for(i=0;i<nc;i++) {
			for(j=0;j<i;j++) {
				Orthogonalize(normX+i,normY+i,normZ+i,normX[j],normY[j],normZ[j]);
			}
			if((normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i])<FFACTOR_ORTH) {
				normX[i] = normY[i] = normZ[i] = 0;
				nconstraint--;
			}
		}

/*
 * 		Find independent line constraints
 */
		nlc = 0;
		for(i=0;i<nc;i++) {
			for(j=0;j<i;j++) {
				Orthogonalize(lineX+i,lineY+i,lineZ+i,lineX[j],lineY[j],lineZ[j]);
			}
			if((lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i])<FFACTOR_ORTH) {
				lineX[i] = lineY[i] = lineZ[i] = 0;
			}
			else {
				nlc++;
			}
		}

/*
 * 		Do not use mobility matrix if node has one line constraint
 */		
		if (nlc == 1) {
			//node->mobMatrixUse = 0;
			// DO SOMETHING HERE...
			//return 0;
		}

/*
 *  	Velocity is simply proportional to total force per unit length
 */		
		for (i = 0; i < 3; i++)
			mobMatrix[i][i] = 1.0;
    
/*
 *  	Orthogonalize with glide plane constraints
 */
		for(i=0;i<nc;i++) {
			if((normX[i]!=0)||(normY[i]!=0)||(normZ[i]!=0)) {
				
				real8 OrthMatrix[3][3], tmpMob[3][3];
				a = normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i];
				
				OrthMatrix[0][0] = 1.0-normX[i]*normX[i]/a;
				OrthMatrix[0][1] = -normX[i]*normY[i]/a;
				OrthMatrix[0][2] = -normX[i]*normZ[i]/a;
				
				OrthMatrix[1][0] = -normY[i]*normX[i]/a;
				OrthMatrix[1][1] = 1.0-normY[i]*normY[i]/a;
				OrthMatrix[1][2] = -normY[i]*normZ[i]/a;
				
				OrthMatrix[2][0] = -normZ[i]*normX[i]/a;
				OrthMatrix[2][1] = -normZ[i]*normY[i]/a;
				OrthMatrix[2][2] = 1.0-normZ[i]*normZ[i]/a;
					  
				Matrix33Mult33(OrthMatrix, mobMatrix, tmpMob);
				
				for (j = 0; j < 3; j++)
					for (k = 0; k < 3; k++)
						mobMatrix[j][k] = tmpMob[j][k];
			}
		}


/* 
 * 		Any dislocation with glide plane not {111} type 
 * 		can only move along its length. This rule includes 
 * 		LC junction which is on {100} plane
 */
#if _ENABLE_LINE_CONSTRAINT
		
		if(nlc==1) {
			
			/* only one line constraint */
			for(i=0;i<nc;i++) {
				if((lineX[i]!=0)||(lineY[i]!=0)||(lineZ[i]!=0)) { 
					lc[0] = lineX[i];
					lc[1] = lineY[i];
					lc[2] = lineZ[i];
					break;
				}
			}
			
			real8 ProjMatrix[3][3], tmpMob[3][3];
			Vec3TransposeAndMult(lc, ProjMatrix);
			
			Matrix33Mult33(ProjMatrix, mobMatrix, tmpMob);
				
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					mobMatrix[j][k] = tmpMob[j][k];

			if (nconstraint<=0) {
				//Fatal("MobilityLaw_FCC_0: nconstraint <= 0, nlc = 1 is impossible!");
				
				/* set velocity to zero if line is not on every plane */
				for (j = 0; j < 3; j++)
					for (k = 0; k < 3; k++)
						mobMatrix[j][k] = 0.0;
				
			} else if(nconstraint>=1) {
				/* a few plane constraints and one line constraint */
				for(i=0;i<nc;i++) {
					normdotlc = normx[i]*lc[0] + normy[i]*lc[1] + normz[i]*lc[2];
					if(fabs(normdotlc)>FFACTOR_ORTH) {
						
						/* set velocity to zero if line is not on every plane */
						for (j = 0; j < 3; j++)
							for (k = 0; k < 3; k++)
								mobMatrix[j][k] = 0.0;
						
						break;
					}
				}
			}
		} else if (nlc>=2) {
			
			/* Velocity is zero when # of independnet lc 
			 * constratins is equal to or more than 2 */
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					mobMatrix[j][k] = 0.0;
		}
#endif
	
		return 0;
}
#endif

#ifdef _PRECOMPUTE_GLIDE
/*-------------------------------------------------------------------------
 *
 *      Function:     Mobility_FCC_0_pre_matrix
 *      Description:  This function precomputes a nodal mobility matrix 
 *                    operator lumping the glide and line constraints
 *                    to evaluate the mobility from forces.
 *                    Note: The B matrix (related to arms length) is 
 *                    approximated as a constant during time-integration
 *                    when using option _PRECOMPUTE_BMATRIX.
 * 
 *-----------------------------------------------------------------------*/
int Mobility_FCC_0_pre_matrix(Home_t *home, Node_t *node)
{
		int numNonZeroLenSegs = 0;
		Param_t *param;
		real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
		int i, j, k, nc, nconstraint, nlc;
		real8 normX[100], normY[100], normZ[100];
		real8 normx[100], normy[100], normz[100];
		real8 lineX[100], lineY[100], lineZ[100];
		real8 a, b;
		real8 dx, dy, dz, lx, ly, lz, lr, LtimesB;
		real8 lcx, lcy, lcz, normdotlc, lc[3];
		Node_t *nbr;
		real8 MobScrew, MobEdge, Mob;
		real8 bx, by, bz, br, dangle;
		real8 nForce[3];

		param = home->param;

		MobScrew = param->MobScrew;
		MobEdge  = param->MobEdge;

		nc = node->numNbrs;
		
/*
 *  	Initialize mobility matrix operator
 */		
		node->mobMatrixUse = 1;
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				node->mobMatrix[i][j] = 0.0;
  
/*
 *  	If node is 'pinned' in place by constraints, or the node has any arms 
 *  	with a burgers vector that has explicitly been set to be sessile (via
 *  	control file inputs), the node may not be moved so just zero the velocity
 *  	and return
 */
		if ((node->constraint == PINNED_NODE) || NodeHasSessileBurg(home, node)) {
			return(0);
		}

/*
 *  	It's possible this function was called for a node which had only zero-
 *  	length segments (during SplitSurfaceNodes() for example).  If that is
 *  	the case, just set the velocity to zero and return.
 */
		for (i = 0; i < nc; i++) {
			if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) continue;
			dx = node->x - nbr->x;
			dy = node->y - nbr->y;
			dz = node->z - nbr->z;
			if ((dx*dx + dy*dy + dz*dz) > 1.0e-12) {
				numNonZeroLenSegs++;
			}
		}

		if (numNonZeroLenSegs == 0) {
			return(0);
		}


/*
 *      Copy glide plane constraints and determine
 *      line constraints.
 */
		for(i=0;i<nc;i++) {
			
			normX[i] = normx[i] = node->nx[i];
			normY[i] = normy[i] = node->ny[i];
			normZ[i] = normz[i] = node->nz[i];

/*
 *      	If needed, rotate the glide plane normals from the
 *      	laboratory frame to the crystal frame.
 */
			if (param->useLabFrame) {
				real8 normTmp[3] = {normX[i], normY[i], normZ[i]};
				real8 normRot[3];

				Matrix33Vector3Multiply(home->rotMatrixInverse, normTmp, normRot);

				normX[i] = normRot[0]; normY[i] = normRot[1]; normZ[i] = normRot[2];
				normx[i] = normRot[0]; normy[i] = normRot[1]; normz[i] = normRot[2];
			}

			if ( (fabs(fabs(normX[i]) - fabs(normY[i])) > FFACTOR_NORMAL) ||
                 (fabs(fabs(normY[i]) - fabs(normZ[i])) > FFACTOR_NORMAL) ) {
				
				/* not {111} plane */
				if ((nbr=GetNeighborNode(home,node,i)) == (Node_t *)NULL) {
					Fatal("Neighbor not found at %s line %d\n",__FILE__,__LINE__);
				}
				
				lineX[i] = nbr->x - node->x;
				lineY[i] = nbr->y - node->y; 
				lineZ[i] = nbr->z - node->z;
				ZImage (param, lineX+i, lineY+i, lineZ+i);

/*
 *          	If needed, rotate the line sense from the laboratory frame to
 *          	the crystal frame.
 */
				if (param->useLabFrame) {
					real8 lDir[3] = {lineX[i], lineY[i], lineZ[i]};
					real8 lDirRot[3];

					Matrix33Vector3Multiply(home->rotMatrixInverse, lDir, lDirRot);

					lineX[i] = lDirRot[0];
					lineY[i] = lDirRot[1];
					lineZ[i] = lDirRot[2];
				}
				
			} else { 
				/* no line constraint */
				lineX[i] = lineY[i] = lineZ[i] = 0;
			}
		}
    
/*
 * 		Normalize glide plane normal vectors 
 *		and lc line vectors
 */
		for(i=0;i<nc;i++) {
			
			a=sqrt(normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i]);
			b=sqrt(lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i]);

			if(a>0) {
				normX[i]/=a;
				normY[i]/=a;
				normZ[i]/=a;

				normx[i]/=a;
				normy[i]/=a;
				normz[i]/=a;
			}
			if(b>0) {
				lineX[i]/=b;
				lineY[i]/=b;
				lineZ[i]/=b;
			}
		}

/*
 * 		Find independent glide constraints
 */
		nconstraint = nc;
		for(i=0;i<nc;i++) {
			for(j=0;j<i;j++) {
				Orthogonalize(normX+i,normY+i,normZ+i,normX[j],normY[j],normZ[j]);
			}
			if((normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i])<FFACTOR_ORTH) {
				normX[i] = normY[i] = normZ[i] = 0;
				nconstraint--;
			}
		}

/*
 * 		Find independent line constraints
 */
		nlc = 0;
		for(i=0;i<nc;i++) {
			for(j=0;j<i;j++) {
				Orthogonalize(lineX+i,lineY+i,lineZ+i,lineX[j],lineY[j],lineZ[j]);
			}
			if((lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i])<FFACTOR_ORTH) {
				lineX[i] = lineY[i] = lineZ[i] = 0;
			}
			else {
				nlc++;
			}
		}

/*
 * 		Do not use mobility matrix if node has line constraints
 */		
		if (nlc > 0) {
			node->mobMatrixUse = 0;
			return(0);
		}

#if _PRECOMPUTE_BMATRIX
/*
 * 		Find total dislocation length times drag coefficent (LtimesB)
 */
		LtimesB=0;
		for(j=0;j<nc;j++) {
			
			if ((nbr=GetNeighborNode(home,node,j)) == (Node_t *)NULL) continue;
			dx=nbr->x - node->x;
			dy=nbr->y - node->y;
			dz=nbr->z - node->z;
			ZImage (param, &dx, &dy, &dz);

/*
 *      	If needed, rotate the line sense from the laboratory frame to
 *      	the crystal frame.
 */
			if (param->useLabFrame) {
				real8 dTmp[3] = {dx, dy, dz};
				real8 dRot[3];

				Matrix33Vector3Multiply(home->rotMatrixInverse, dTmp, dRot);

				dx = dRot[0]; dy = dRot[1]; dz = dRot[2];
			}

			lr=sqrt(dx*dx+dy*dy+dz*dz);
        
			if (lr==0) {
				/* zero arm segment can happen after node split 
				 * it is OK to have plane normal vector == 0
				 * Skip (do nothing)
				 */
			} else {
			
				if((node->nx[j]==0)&&(node->ny[j]==0)&&(node->nz[j]==0)) {
				   
					printf("Mobility_FCC_0: (%d,%d) glide plane norm = 0\n"
					       "for segment with nonzero length lr = %e!\n",
					       node->myTag.domainID, node->myTag.index, lr);
				}

				lx=dx/lr; ly=dy/lr; lz=dz/lr;

				bx = node->burgX[j];
				by = node->burgY[j];
				bz = node->burgZ[j];
/*
 * 				If needed, rotate the burgers vector from the laboratory frame to
 *         		the crystal frame.
 */
				if (param->useLabFrame) {
					real8 bTmp[3] = {bx, by, bz};
					real8 bRot[3];

					Matrix33Vector3Multiply(home->rotMatrixInverse, bTmp, bRot);

					bx = bRot[0]; by = bRot[1]; bz = bRot[2];
				}

				br = sqrt(bx*bx+by*by+bz*bz);
				bx/=br; by/=br; bz/=br; /* unit vector along Burgers vector */

				dangle = fabs(bx*lx+by*ly+bz*lz);
				Mob=MobEdge+(MobScrew-MobEdge)*dangle;

				LtimesB+=(lr/Mob);
			}
		}
		LtimesB/=2;

/*
 *  	Velocity is simply proportional to total force per unit length
 */	
		for (i = 0; i < 3; i++)
			node->mobMatrix[i][i] = 1.0/LtimesB;

#else
/*
 *  	Velocity is simply proportional to total force per unit length
 */		
		for (i = 0; i < 3; i++)
			node->mobMatrix[i][i] = 1.0;
#endif
    
/*
 *  	Orthogonalize with glide plane constraints
 */
		for(i=0;i<nc;i++) {
			if((normX[i]!=0)||(normY[i]!=0)||(normZ[i]!=0)) {
				
				real8 OrthMatrix[3][3], tmpMob[3][3];
				a = normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i];
				
				OrthMatrix[0][0] = 1.0-normX[i]*normX[i]/a;
				OrthMatrix[0][1] = -normX[i]*normY[i]/a;
				OrthMatrix[0][2] = -normX[i]*normZ[i]/a;
				
				OrthMatrix[1][0] = -normY[i]*normX[i]/a;
				OrthMatrix[1][1] = 1.0-normY[i]*normY[i]/a;
				OrthMatrix[1][2] = -normY[i]*normZ[i]/a;
				
				OrthMatrix[2][0] = -normZ[i]*normX[i]/a;
				OrthMatrix[2][1] = -normZ[i]*normY[i]/a;
				OrthMatrix[2][2] = 1.0-normZ[i]*normZ[i]/a;
					  
				Matrix33Mult33(OrthMatrix, node->mobMatrix, tmpMob);
				
				for (j = 0; j < 3; j++)
					for (k = 0; k < 3; k++)
						node->mobMatrix[j][k] = tmpMob[j][k];
			}
		}


/* 
 * 		Any dislocation with glide plane not {111} type 
 * 		can only move along its length. This rule includes 
 * 		LC junction which is on {100} plane
 */
#if _ENABLE_LINE_CONSTRAINT
		
		if(nlc==1) {
			
			/* only one line constraint */
			for(i=0;i<nc;i++) {
				if((lineX[i]!=0)||(lineY[i]!=0)||(lineZ[i]!=0)) { 
					lc[0] = lineX[i];
					lc[1] = lineY[i];
					lc[2] = lineZ[i];
					break;
				}
			}
			
			real8 ProjMatrix[3][3], tmpMob[3][3];
			Vec3TransposeAndMult(lc, ProjMatrix);
			
			Matrix33Mult33(ProjMatrix, node->mobMatrix, tmpMob);
				
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					node->mobMatrix[j][k] = tmpMob[j][k];

			if (nconstraint<=0) {	
				Fatal("MobilityLaw_FCC_0: nconstraint <= 0, nlc = 1 is impossible!");
			} else if(nconstraint>=1) {
				/* a few plane constraints and one line constraint */
				for(i=0;i<nc;i++) {
					normdotlc = normx[i]*lc[0] + normy[i]*lc[1] + normz[i]*lc[2];
					if(fabs(normdotlc)>FFACTOR_ORTH) {
						
						/* set velocity to zero if line is not on every plane */
						for (j = 0; j < 3; j++)
							for (k = 0; k < 3; k++)
								node->mobMatrix[j][k] = 0.0;
						
						break;
					}
				}
			}
		} else if (nlc>=2) {
			
			/* Velocity is zero when # of independnet lc 
			 * constratins is equal to or more than 2 */
			for (j = 0; j < 3; j++)
				for (k = 0; k < 3; k++)
					node->mobMatrix[j][k] = 0.0;
		}
#endif
	
		return(0);
}

/*-------------------------------------------------------------------------
 *
 *      Function:     Mobility_FCC_0_matrix
 *      Description:  This function computes the velocity of a node from
 *                    the force in the case of FCC crystals after the
 *                    glide and line constraints associated with the node
 *                    have been lumped into the nodal mobility matrix 
 *                    operator using function Mobility_FCC_0_pre_matrix
 * 
 *-----------------------------------------------------------------------*/
int Mobility_FCC_0_matrix(Home_t *home, Node_t *node)
{
		int   j, k, nc, numNonZeroLenSegs;
		real8 nForce[3], nVel[3];
		real8 VelxNode, VelyNode, VelzNode;
		real8 dx, dy, dz, lx, ly, lz, lr, LtimesB;
		real8 MobScrew, MobEdge, Mob;
		real8 bx, by, bz, br, dangle;
		Node_t  *nbr;
		Param_t *param;
		
		param = home->param;
		
		nForce[0] = node->fX;
		nForce[1] = node->fY;
		nForce[2] = node->fZ;
		
/*
 *  	If needed, rotate the force vector from the laboratory frame to the
 *  	crystal frame
 */
		if (param->useLabFrame) {
			real8 rotForce[3];
			Matrix33Vector3Multiply(home->rotMatrixInverse, nForce, rotForce);
			VECTOR_COPY(nForce, rotForce);
		}

#if !_PRECOMPUTE_BMATRIX
/*
 * 		Find total dislocation length times drag coefficent (LtimesB)
 */	
		MobScrew = param->MobScrew;
		MobEdge  = param->MobEdge;
		
		nc = node->numNbrs;
		numNonZeroLenSegs = 0;
		LtimesB = 0;
		
		for(j=0;j<nc;j++) {
			
			if ((nbr=GetNeighborNode(home,node,j)) == (Node_t *)NULL) continue;
			dx=nbr->x - node->x;
			dy=nbr->y - node->y;
			dz=nbr->z - node->z;
			ZImage (param, &dx, &dy, &dz);

/*
 *      	If needed, rotate the line sense from the laboratory frame to
 *      	the crystal frame.
 */
			if (param->useLabFrame) {
				real8 dTmp[3] = {dx, dy, dz};
				real8 dRot[3];

				Matrix33Vector3Multiply(home->rotMatrixInverse, dTmp, dRot);

				dx = dRot[0]; dy = dRot[1]; dz = dRot[2];
			}

			lr=sqrt(dx*dx+dy*dy+dz*dz);
        
			if (lr > 1.e-12) numNonZeroLenSegs++;
				
			if (lr != 0) {
				
				if((node->nx[j]==0)&&(node->ny[j]==0)&&(node->nz[j]==0)) {
				   
					printf("Mobility_FCC_0: (%d,%d) glide plane norm = 0\n"
					       "for segment with nonzero length lr = %e!\n",
					       node->myTag.domainID, node->myTag.index, lr);
				}

				lx=dx/lr; ly=dy/lr; lz=dz/lr;

				bx = node->burgX[j];
				by = node->burgY[j];
				bz = node->burgZ[j];
/*
 * 				If needed, rotate the burgers vector from the laboratory frame to
 *         		the crystal frame.
 */
				if (param->useLabFrame) {
					real8 bTmp[3] = {bx, by, bz};
					real8 bRot[3];

					Matrix33Vector3Multiply(home->rotMatrixInverse, bTmp, bRot);

					bx = bRot[0]; by = bRot[1]; bz = bRot[2];
				}

				br = sqrt(bx*bx+by*by+bz*bz);
				bx/=br; by/=br; bz/=br; /* unit vector along Burgers vector */

				dangle = fabs(bx*lx+by*ly+bz*lz);
				Mob=MobEdge+(MobScrew-MobEdge)*dangle;

				LtimesB+=(lr/Mob);
			}
		}
		LtimesB/=2;
		
		nForce[0] /= LtimesB;
		nForce[1] /= LtimesB;
        nForce[2] /= LtimesB;
          
/*
 *  	It's possible this function was called for a node which had only zero-
 *  	length segments (during SplitSurfaceNodes() for example).  If that is
 *  	the case, just set the velocity to zero and return.
 */		
		if (numNonZeroLenSegs == 0) {
			node->vX = 0.0;
			node->vY = 0.0;
			node->vZ = 0.0;
			return(0);
		}
#endif
        
/*
 *  	Calculate the velocity
 */
		
		Matrix33Vector3Multiply(node->mobMatrix, nForce, nVel);
		VelxNode = nVel[0];
		VelyNode = nVel[1];
		VelzNode = nVel[2];
		
/*
 *  	If needed, rotate the velocity vector back to the laboratory frame
 *  	from the crystal frame
 */
		if (param->useLabFrame) {
			real8 vTmp[3] = {VelxNode, VelyNode, VelzNode};
			real8 vRot[3];

			Matrix33Vector3Multiply(home->rotMatrix, vTmp, vRot);

			VelxNode = vRot[0];
			VelyNode = vRot[1];
			VelzNode = vRot[2];
		}

		node->vX = VelxNode;
		node->vY = VelyNode;
		node->vZ = VelzNode;
		
		return(0);
}
#endif

/*-------------------------------------------------------------------------
 *
 *      Function:     Mobility_FCC_0_original
 *      Description:  This function computes the velocity of a node from
 *                    the force in the case of FCC crystals.
 *
 *-----------------------------------------------------------------------*/
int Mobility_FCC_0_original(Home_t *home, Node_t *node)
{
    int numNonZeroLenSegs = 0;
    Param_t *param;
    real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
    int i, j, nc, nconstraint, nlc;
    real8 normX[100], normY[100], normZ[100], normx[100], normy[100], normz[100];
    real8 lineX[100], lineY[100], lineZ[100];
    real8 a, b;
    real8 dx, dy, dz, lx, ly, lz, lr, LtimesB;
    real8 lcx, lcy, lcz, normdotlc;
    Node_t *nbr;
    real8 MobScrew, MobEdge, Mob;
    real8 bx, by, bz, br, dangle;
    real8 nForce[3];

    param = home->param;

    MobScrew = param->MobScrew;
    MobEdge  = param->MobEdge;
    
    nc = node->numNbrs;
    
/*
 *  If node is 'pinned' in place by constraints, or the node has any arms 
 *  with a burgers vector that has explicitly been set to be sessile (via
 *  control file inputs), the node may not be moved so just zero the velocity
 *  and return
 */
    if ((node->constraint == PINNED_NODE) ||
        NodeHasSessileBurg(home, node))
    {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }

/*
 *  It's possible this function was called for a node which had only zero-
 *  length segments (during SplitSurfaceNodes() for example).  If that is
 *  the case, just set the velocity to zero and return.
 */
    for (i = 0; i < nc; i++) {
        if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) continue;
        dx = node->x - nbr->x;
        dy = node->y - nbr->y;
        dz = node->z - nbr->z;
        if ((dx*dx + dy*dy + dz*dz) > 1.0e-12) {
            numNonZeroLenSegs++;
        }
    }

    if (numNonZeroLenSegs == 0) {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }


    /* copy glide plane constraints and determine line constraints */
    for(i=0;i<nc;i++)
    {
        normX[i] = normx[i] = node->nx[i];
        normY[i] = normy[i] = node->ny[i];
        normZ[i] = normz[i] = node->nz[i];

/*
 *      If needed, rotate the glide plane normals from the
 *      laboratory frame to the crystal frame.
 */
        if (param->useLabFrame) {
            real8 normTmp[3] = {normX[i], normY[i], normZ[i]};
            real8 normRot[3];

            Matrix33Vector3Multiply(home->rotMatrixInverse, normTmp, normRot);

            normX[i] = normRot[0]; normY[i] = normRot[1]; normZ[i] = normRot[2];
            normx[i] = normRot[0]; normy[i] = normRot[1]; normz[i] = normRot[2];
        }

        if ( (fabs(fabs(normX[i]) - fabs(normY[i])) > FFACTOR_NORMAL) ||
             (fabs(fabs(normY[i]) - fabs(normZ[i])) > FFACTOR_NORMAL) )
        { /* not {111} plane */
            if ((nbr=GetNeighborNode(home,node,i)) == (Node_t *)NULL) {
                Fatal("Neighbor not found at %s line %d\n",__FILE__,__LINE__);
            }
            lineX[i] = nbr->x - node->x;
            lineY[i] = nbr->y - node->y; 
            lineZ[i] = nbr->z - node->z;
            ZImage (param, lineX+i, lineY+i, lineZ+i);

/*
 *          If needed, rotate the line sense from the laboratory frame to
 *          the crystal frame.
 */
            if (param->useLabFrame) {
                real8 lDir[3] = {lineX[i], lineY[i], lineZ[i]};
                real8 lDirRot[3];

                Matrix33Vector3Multiply(home->rotMatrixInverse, lDir, lDirRot);

                lineX[i] = lDirRot[0];
                lineY[i] = lDirRot[1];
                lineZ[i] = lDirRot[2];
            }
	}
	else
	{ /* no line constraint */
	    lineX[i] = lineY[i] = lineZ[i] = 0;
	}
    }
    
    /* normalize glide plane normal vectors and lc line vectors*/
    for(i=0;i<nc;i++)
    {
        a=sqrt(normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i]);
	b=sqrt(lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i]);

        if(a>0)
        {
            normX[i]/=a;
            normY[i]/=a;
            normZ[i]/=a;

            normx[i]/=a;
            normy[i]/=a;
            normz[i]/=a;
        }
        if(b>0)
        {
            lineX[i]/=b;
            lineY[i]/=b;
            lineZ[i]/=b;
        }
    }


    /* Find independent glide constraints */ 
    nconstraint = nc;
    for(i=0;i<nc;i++)
    {
        for(j=0;j<i;j++)
        {
            Orthogonalize(normX+i,normY+i,normZ+i,normX[j],normY[j],normZ[j]);
        }
        if((normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i])<FFACTOR_ORTH)
        {
            normX[i] = normY[i] = normZ[i] = 0;
            nconstraint--;
        }
    }

    /* Find independent line constraints */
    nlc = 0;
    for(i=0;i<nc;i++)
    {
        for(j=0;j<i;j++)
        {
            Orthogonalize(lineX+i,lineY+i,lineZ+i,lineX[j],lineY[j],lineZ[j]);
        }
        if((lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i])<FFACTOR_ORTH)
        {
            lineX[i] = lineY[i] = lineZ[i] = 0;
        }
        else
        {
            nlc ++;
        }
    }

    /* find total dislocation length times drag coefficent (LtimesB)*/
    LtimesB=0;
    for(j=0;j<node->numNbrs;j++)
    {
        if ((nbr=GetNeighborNode(home,node,j)) == (Node_t *)NULL) continue;
        dx=nbr->x - node->x;
        dy=nbr->y - node->y;
        dz=nbr->z - node->z;
        ZImage (param, &dx, &dy, &dz) ;

/*
 *      If needed, rotate the line sense from the laboratory frame to
 *      the crystal frame.
 */
        if (param->useLabFrame) {
            real8 dTmp[3] = {dx, dy, dz};
            real8 dRot[3];

            Matrix33Vector3Multiply(home->rotMatrixInverse, dTmp, dRot);

            dx = dRot[0]; dy = dRot[1]; dz = dRot[2];
        }

        lr=sqrt(dx*dx+dy*dy+dz*dz);
        
        if (lr==0)
        { /* zero arm segment can happen after node split 
           * it is OK to have plane normal vector == 0
           * Skip (do nothing)
           */
        }
        else 
        {
           if((node->nx[j]==0)&&(node->ny[j]==0)&&(node->nz[j]==0))
           {
              printf("Mobility_FCC_0: (%d,%d) glide plane norm = 0\n"
                     "for segment with nonzero length lr = %e!\n",
                     node->myTag.domainID, node->myTag.index, lr);
           }

           lx=dx/lr; ly=dy/lr; lz=dz/lr;

           bx = node->burgX[j];
           by = node->burgY[j];
           bz = node->burgZ[j];
/*
 *         If needed, rotate the burgers vector from the laboratory frame to
 *         the crystal frame.
 */
           if (param->useLabFrame) {
               real8 bTmp[3] = {bx, by, bz};
               real8 bRot[3];

               Matrix33Vector3Multiply(home->rotMatrixInverse, bTmp, bRot);

               bx = bRot[0]; by = bRot[1]; bz = bRot[2];
           }

           br = sqrt(bx*bx+by*by+bz*bz);
           bx/=br; by/=br; bz/=br; /* unit vector along Burgers vector */

           dangle = fabs(bx*lx+by*ly+bz*lz);
           Mob=MobEdge+(MobScrew-MobEdge)*dangle;

           LtimesB+=(lr/Mob);
	}
    }
    LtimesB/=2;

    nForce[0] = node->fX;
    nForce[1] = node->fY;
    nForce[2] = node->fZ;

/*
 *  If needed, rotate the force vector from the laboratory frame to the
 *  crystal frame
 */
    if (param->useLabFrame) {
        real8 rotForce[3];
        Matrix33Vector3Multiply(home->rotMatrixInverse, nForce, rotForce);
        VECTOR_COPY(nForce, rotForce);
    }

    /* Velocity is simply proportional to total force per unit length */
    VelxNode = nForce[0]/LtimesB;
    VelyNode = nForce[1]/LtimesB;
    VelzNode = nForce[2]/LtimesB;
    

    /* Orthogonalize with glide plane constraints */
    if (nconstraint>=3){
		VelxNode = VelyNode = VelzNode = 0.0;
	} else {
		for(i=0;i<nc;i++)
		{
			if((normX[i]!=0)||(normY[i]!=0)||(normZ[i]!=0))
			{
			 Orthogonalize(&VelxNode,&VelyNode,&VelzNode,
							 normX[i],normY[i],normZ[i]);
			}
		}
	}
	


    /* Any dislocation with glide plane not {111} type can only move along its length
     * This rule includes LC junction which is on {100} plane
     */
#if _ENABLE_LINE_CONSTRAINT
    if(nlc==1)
    { /* only one line constraint */
        for(i=0;i<nc;i++)
        {
            if((lineX[i]!=0)||(lineY[i]!=0)||(lineZ[i]!=0))
            { 
   	        lcx = lineX[i];
	        lcy = lineY[i];
	        lcz = lineZ[i];
                break;
            }
        }

        /* project velocity along line */
        Veldotlcr = VelxNode*lcx+VelyNode*lcy+VelzNode*lcz;
        VelxNode = Veldotlcr*lcx;
        VelyNode = Veldotlcr*lcy;
        VelzNode = Veldotlcr*lcz;

	if (nconstraint<=0)
	{	
            Fatal("MobilityLaw_FCC_0: nconstraint <= 0, nlc = 1 is impossible!");
        }
        else if(nconstraint>=1)
  	{ /* a few plane constraints and one line constraint */
            for(i=0;i<nc;i++)
            {
		normdotlc = normx[i]*lcx + normy[i]*lcy + normz[i]*lcz;
		if(fabs(normdotlc)>FFACTOR_ORTH)
		{
                    /* set velocity to zero if line is not on every plane */
                    VelxNode = VelyNode = VelzNode = 0;
		    break;
		}
                else
                {
                   /* do nothing. Skip */
                }
	    }
	}
    }
    else if (nlc>=2)
    {
        /* Velocity is zero when # of independnet lc constratins is equal to or more than 2 */
        VelxNode = VelyNode = VelzNode = 0;
    }
    else
    { 
        /* nlc == 0, do nothing */
    }
#endif

/*
 *  If needed, rotate the velocity vector back to the laboratory frame
 *  from the crystal frame
 */
    if (param->useLabFrame) {
        real8 vTmp[3] = {VelxNode, VelyNode, VelzNode};
        real8 vRot[3];

        Matrix33Vector3Multiply(home->rotMatrix, vTmp, vRot);

        VelxNode = vRot[0];
        VelyNode = vRot[1];
        VelzNode = vRot[2];
    }

    node->vX = VelxNode;
    node->vY = VelyNode;
    node->vZ = VelzNode;

    return(0);
}
