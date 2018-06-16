/**************************************************************************
 *
 *      Module:      ImplicitIntegratorSub0.c
 *      Description: Used with the SubcycleIntegratorNodeB function (node based subcycling).
 *                   Implements a numerical integrator using
 *                   the implicitly-solved trapezoid integration method.
 *					 This couples the trapezoid method with the Newton-Raphson
 *					 method for nonlinear system solution. The difference
 *					 between this function and ImplicitIntegrator is that 
 *					 it only time integrates the group 0 nodes. 
 *
 ***************************************************************************/
#ifdef _IMPLICIT

#include "Home.h"
#include "sys/stat.h"
#include "sys/types.h"
#include "Implicit.h"
#include "cs.h"

#ifdef _LIS
#include "lis.h"
#endif

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

#define DIR_TIMESTEP_ERROR "timestep_error"
#define DIR_NODAL_TIMESTEP "timestep"

/*------------------------------------------------------------------------
 *
 *      Function:    PreserveNodalDataSub0
 *      Description: Both old and new values for certain nodal
 *                   data items are required during timestep
 *                   integration and calculating plastic strain.
 *                   This function copies the values for specified
 *                   items into the appropriate variables.
 *                   
 *-----------------------------------------------------------------------*/
static void PreserveNodalDataSub0(Home_t *home)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if (node->subgroup!=0) continue;
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
 *      Function:    FlagGroup0Nodes
 *      Description: Loops over all nodes and toggles the NODE_RESET_FORCES
 *					 for all group 0 nodes.
 *
 *-----------------------------------------------------------------------*/

void FlagGroup0Nodes(Home_t *home) 
{
		int		i;
		Node_t	*node;
		for (i=0; i < home->newNodeKeyPtr; i++) {
	        node = home->nodeKeys[i];
	        if (node == (Node_t *)NULL) continue;
			node->flags &= (~NODE_RESET_FORCES);
			if (node->subgroup==0) {
				node->flags |= NODE_RESET_FORCES;
			}
		}
}


/*------------------------------------------------------------------------
 *
 *      Function:    ImplicitIntegratorSub0
 *      Description: Implements a numerical timestep integrator using
 *                   the implicitly-solved Trapezoid integration method.
 *					 Only integrates group 0 nodes.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
static int Implicit_Step(Home_t *home , real8 newDT)
{
        int     i, j, convergent, maxIterations, tmp;
        int 	ind1, ind2, ind3;
        int     iter, globalIterError, mobIterError;
        int     dumpErrorData = 1, doAll = 1;
		int 	changedJ, errInc, count, NJ, Ndof;
		real8   errMax, globalErrMax, oldglobalErrMax, drn, relerrMax;
        real8   oldx, oldy, oldz, globalRelErrMax;
        real8   localVals[3], globalVals[3];
		real8 	Q[3][3], Qinv[3][3];
		real8	rtmp[3], rout[3], vtmp[3], vout[3];
		real8	thisx, thisy, thisz;
		real8 	errdX, errdY, errdZ, f1, f2, f3;
		real8	dx, dy, dz, xtmp, ytmp, ztmp;
        Node_t  *node;
        Param_t *param;

		param = home->param;
		int N = home->newNodeKeyPtr;

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
 * 		Calculate coordinate systems and transformation matrices for every
 * 		node.
 */		
		for (i=0; i < N; i++) {
	        node = home->nodeKeys[i];
	        if (node == (Node_t *)NULL) continue;
			node->flags |= UPDATE_NODE_J;
		}
		Ndof = 0;
		AssignCoordSys(home, &Ndof, 0);
		
/*
 *		Set the initial guess to the current nodal positions
 */
		for (i=0; i < N; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if (node->ndof==0) continue;
			
			node->x = node->oldx;
			node->y = node->oldy;
			node->z = node->oldz;
		}
		
		mobIterError  = 0;
		errInc        = 0;
		globalVals[0] = 0.0;
		globalVals[1] = 0.0;
		globalVals[2] = 0.0;
		
/*
 *		Calculate nodal forces and velocities
 */
		FlagGroup0Nodes(home);
    	NodeForce(home, PARTIAL);
        mobIterError = CalcNodeVelocities(home, 0, !doAll);
        CommSendVelocity(home);
		if (mobIterError) return(0);
		
/*
 *      Loop until we converge on a time step.  
 */
        convergent    = 0;
        maxIterations = 20;
		
/*
 *		Calculate the Jacobian
 */ 
		Jmobtrp = cs_spalloc(Ndof, Ndof, 2*Ndof, 1, 1);
		UpdateJacobian(home,1,&changedJ,Jmobtrp);  // Jmobtrp gets updated here. 1 means update all nodes. changedJ is useless.

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
				FlagGroup0Nodes(home);
				NodeForce(home, PARTIAL);
				mobIterError = CalcNodeVelocities(home, 0, !doAll);
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
				if ( node->ndof==0) continue;

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
				
				if (thiserr > param->rTolth/2e1) {
					if (drn > param->rTolth / param->rTolrel) {
						relerrMax = MAX( relerrMax , thiserr/drn);
					} else {
						relerrMax = 2.0*param->rTolrel;
					}
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
 *			If the error is within the tolerance, we've reached
 *			convergence so we can accept this deltaT.  Otherwise
 *			continue iterating, unless max number of iterations has reached.
 *			If it has, exit the for loop and cut the time step.  
 *			Note: we need to reposition both local nodes and ghost nodes!
 */
			if (globalErrMax < param->rTol/2e1 && globalRelErrMax < param->rTolrel/2e1
			                               && iter > 0) {
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
 *      Function:    ImplicitIntegratorSub0
 *      Description: Implements a numerical timestep integrator using
 *                   the implicitly-solved Trapezoid integration method.
 *					 Only integrates group 0 nodes.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
void ImplicitIntegratorSub0(Home_t *home)
{
		int     i, convergent, incrDelta;
		int     globalIterError, mobIterError;
		int     doAll = 1, N = home->newNodeKeyPtr;
		real8   errMax, localVals[2], globalVals[2];
		real8   newDT , oldx  , oldy  , oldz     , drn;
		real8	errx  , erry  , errz  , thiserr  , factor;
		real8   errdX , errdY , errdZ , relerrMax, globalErrMax;
		real8   tmp1  , tmp2  , tmp3  , tmp4     , globalRelErrMax;
		Node_t  *node;
		Param_t *param;

		param = home->param;
		
		int    *col = malloc(N * sizeof(int     ));
		real8  **X0 = malloc(N * sizeof(real8 * ));
		real8  **X1 = malloc(N * sizeof(real8 * ));
		real8  **X2 = malloc(N * sizeof(real8 * ));
		
		for ( i = 0 ; i < N ; i++ ) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			X0[i] = malloc( 3 * sizeof(real8 ));
			X1[i] = malloc( 3 * sizeof(real8 ));
			X2[i] = malloc( 3 * sizeof(real8 ));
			
			X0[i][0] = node->x;
			X0[i][1] = node->y;
			X0[i][2] = node->z;
		}

/*
 * 		Set the initial time step.
 */
		newDT = MIN(param->maxDT, param->nextDT);
		if (newDT <= 0.0) newDT = param->maxDT;
		param->deltaTT = newDT;

/*
 *      Preserve certain items of nodal data for later use.
 */
        PreserveNodalDataSub0(home);

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
		
/*
 *			Taking one step of size dt and storing the positions in X1
 */
			convergent = Implicit_Step(home , newDT);
			if (!convergent) {
				iTry1++;
			/*	if (iTry1==1) {
					int ioGroup = home->ioGroupNum;
					char debug_data [50];
					sprintf (debug_data, "first%d", home->cycle);
					WriteAtomEye(home, debug_data, ioGroup, 1, 1, 1);
				}*/
			
				cutTimeStepSize(home , &newDT , &incrDelta);
				continue;
			}
		//	continue;
			
			iTry2++;
			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				X1[i][0] = node->x;
				X1[i][1] = node->y;
				X1[i][2] = node->z;
			}
			
/*
 *			Taking two RK steps of size dt/2
 *			First take one step of size dt/2 starting form the old positions
 */
			convergent = Implicit_Step(home , newDT*0.5);
			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				node->oldx = node->x;   node->currvX = node->vX;
				node->oldy = node->y;   node->currvY = node->vY;
				node->oldz = node->z;   node->currvZ = node->vZ;
			}
			
/*
 *			Now take the second half step and evaluate the truncation error
 */
			convergent = Implicit_Step(home , newDT*0.5);
			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				node->oldx = X0[i][0];  X2[i][0] = node->x;
				node->oldy = X0[i][1];  X2[i][1] = node->y;
				node->oldz = X0[i][2];  X2[i][2] = node->z;
				
				node->currvX = node->oldvX;
				node->currvY = node->oldvY;
				node->currvZ = node->oldvZ;
				
/*
 *				Calculate the truncation error on this node
 */
				errx = X2[i][0] - X1[i][0];
				erry = X2[i][1] - X1[i][1];
				errz = X2[i][2] - X1[i][2];
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
				
			/*	if (iTry2==1) {
					if ((thiserr < param->rTol) && (thiserr < param->rTolth || thiserr/drn < param->rTolrel)) {
						col[i] = 0;
					} else {
						col[i] = 1;
					}
				}*/
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
				cutTimeStepSize(home , &newDT , &incrDelta);
			}
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
		param->deltaTT   = newDT;
		param->realdt    = newDT;
		param->timeStart = param->timeNow;
		
		if (incrDelta) {
            if (param->dtVariableAdjustment) {
                tmp1 = pow(param->dtIncrementFact, param->dtExponent);
                tmp2 = globalErrMax/param->rTol;
                tmp3 = 1.0 / param->dtExponent;
                tmp4 = pow(1.0/(1.0+(tmp1-1.0)*tmp2), tmp3);
                factor = param->dtIncrementFact * tmp4;
				
                param->nextDT = MIN(param->maxDT, newDT*factor);
            } else {
				param->nextDT = MIN(param->maxDT, newDT*param->dtIncrementFact);
            }
        } else {
            param->nextDT = newDT;
        }
		
/*
 *		Freeing memory
 */ 
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			free( X0[i] );
			free( X1[i] );
			free( X2[i] );
		}
		free( X0 );
		free( X1 );
		free( X2 );
		free( col);
		
		return;
		
}

#endif //_IMPLICIT
