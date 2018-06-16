
/**************************************************************************
 *
 *      Module:      TrapezoidIntegratorSub1.c
 *      Description: Used with the Subcycle Integrator function.
 *					 Implements a numerical timestep integrator using
 *                   the Trapezoid integration method.
 *
 ***************************************************************************/
#include "Home.h"
#include "Comm.h"
#include "sys/stat.h"
#include "sys/types.h"

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

void UnpackQMatrices(Node_t *node, real8 Q[3][3], real8 Qinv[3][3]);
void AssignCoordSys(Home_t *home, int *Ndof, int group);
void NodeForceList(Home_t *home, Node_t ***segseglist, int segseglistsize, int reqType );
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
static void PreserveNodalDataSub1(Home_t *home, int items)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
 
			if (node->subgroup!=1) {
				continue;
			}
/*
 *          Save previous nodal position.
 */
            if (items & NODE_POSITION) {
                node->oldx = node->x;
                node->oldy = node->y;
                node->oldz = node->z;
            }

/*
 *          Make a copy of the current velocity without wiping out the
 *          copy of the previous velocity.
 */
            if (items & NODE_CURR_VEL) {
                node->currvX = node->vX;
                node->currvY = node->vY;
                node->currvZ = node->vZ;
            }

/*
 *          Make the current velocity the previous velocity (done at the
 *          end of the timestep integrator)
 */
            if (items & NODE_OLD_VEL) {
                node->oldvX = node->currvX;
                node->oldvY = node->currvY;
                node->oldvZ = node->currvZ;
            }
        }

        node = home->ghostNodeQ;

        while (node != (Node_t *)NULL) {

            if (items & NODE_POSITION) {
                node->oldx = node->x;
                node->oldy = node->y;
                node->oldz = node->z;
            }

            if (items & NODE_CURR_VEL) {
                node->currvX = node->vX;
                node->currvY = node->vY;
                node->currvZ = node->vZ;
            }

            if (items & NODE_OLD_VEL) {
                node->oldvX = node->currvX;
                node->oldvY = node->currvY;
                node->oldvZ = node->currvZ;
            }

            node = node->next;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *	Function:	AdvanceAllNodesSub1
 *	Description:	Reposition all nodes (local and ghost) based
 *                      on the old/current nodal velocities and
 *                      time deltas.
 *
 *			This function assumes current nodal positions and
 *                      velocities have already been preserved in
 *                      the old* variables by a call to PreserveNodalData().
 *
 *      Given the following:
 *
 *      currDT = desired delta time this timestep
 *      oldDT  = delta time used in the previous timestep
 *      currV  = velocity on entry to timestep integrator
 *      oldV   = velocity on entry to timestep integrator on previous step
 *      currP  = current nodal position
 *      newP   = initial guess at position node will end up in after 
 *               this timestep.
 *
 *      The new positions are calculates as:
 *
 *          newP = currP + 0.5 * (currV + oldV) * currDT;
 *
 *-------------------------------------------------------------------------*/
static void AdvanceAllNodesSub1(Home_t *home, real8 oldDT)
{
		int	i;
		real8	x, y, z;
        real8   currDT;
		real8	vtmp[3], currvtmp[3], vout[3], currvout[3], drloc[3], dr[3];
		real8 	Q[3][3], Qinv[3][3];
		Node_t	*node;
		Param_t	*param;

	param = home->param;
	param->realdtsub = param->deltaTTsub;

	currDT = param->realdtsub;

        if (oldDT <= 0.0) oldDT = currDT;

		for (i = 0; i < home->newNodeKeyPtr; i++) {

			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

/*
 *			If this node does not belong to group 1, don't touch it!.
 */
			if (node->subgroup!=1) continue;

/*
 *          If we don't have a value for the previous velocity, assume
 *          previous velocity was same as current velocity.
 */
            if ((node->oldvX == 0.0) && (node->oldvY == 0.0) &&
                (node->oldvZ == 0.0)) {
                node->oldvX = node->currvX;
                node->oldvY = node->currvY;
                node->oldvZ = node->currvZ;
            }

/*
 *				Transform the velocities into local coordinates.
 */
                UnpackQMatrices(node, Q, Qinv);
				vtmp[0] = node->oldvX;
				vtmp[1] = node->oldvY;
				vtmp[2] = node->oldvZ;					
				currvtmp[0] = node->currvX;
				currvtmp[1] = node->currvY;
				currvtmp[2] = node->currvZ;

				Matrix33Vector3Multiply(Qinv, vtmp, vout);
				Matrix33Vector3Multiply(Qinv, currvtmp, currvout);

/*
 *				Calculate the updates as appropriate for the local
 *				coordinates that correspond to degrees of freedom.
 */
				if (node->ndof>0) {
					drloc[0] = (vout[0]+currvout[0])*0.5*currDT;
				} else {
					drloc[0] = 0;
				}
				if (node->ndof>1) {
					drloc[1] = (vout[1]+currvout[1])*0.5*currDT;
				} else {
					drloc[1] = 0;
				}
				if (node->ndof>2) {
					drloc[2] = (vout[2]+currvout[2])*0.5*currDT;
				} else {
					drloc[2] = 0;
				}

/*						
 *				Transform the updates into global coordinates and 
 *				modify the position.
 */
				Matrix33Vector3Multiply(Q, drloc, dr);

                x = node->oldx + dr[0];
                y = node->oldy + dr[1];
                z = node->oldz + dr[2];

				FoldBox(param,&x,&y,&z);
        
				node->x = x;
				node->y = y;
				node->z = z;
	}

/*
 *	Also need to move Ghost nodes
 */
	node = home->ghostNodeQ;

	while (node) {

                if ((node->oldvX == 0.0) && (node->oldvY == 0.0) &&
                    (node->oldvZ == 0.0)) {
                    node->oldvX = node->currvX;
                    node->oldvY = node->currvY;
                    node->oldvZ = node->currvZ;
                }

                x = node->oldx + 0.5 * (node->currvX + node->oldvX) * currDT;
                y = node->oldy + 0.5 * (node->currvY + node->oldvY) * currDT;
                z = node->oldz + 0.5 * (node->currvZ + node->oldvZ) * currDT;

		FoldBox(param,&x,&y,&z);
        
		node->x = x;
		node->y = y;
		node->z = z;

		node = node->next;
	}

	return;    
}

/*------------------------------------------------------------------------
 *
 *      Function:    FlagGroup1Nodes
 *      Description: Loops over all nodes and toggles the NODE_RESET_FORCES
 *					 for all group 1 nodes.
 *
 *-----------------------------------------------------------------------*/

void FlagGroup1Nodes(Home_t *home) 
{
		int		i;
		Node_t	*node;
		for (i=0; i < home->newNodeKeyPtr; i++) {
	        node = home->nodeKeys[i];
	        if (node == (Node_t *)NULL) continue;
			node->flags &= (~NODE_RESET_FORCES);
			if (node->subgroup==1) {
				node->flags |= NODE_RESET_FORCES;
			}
		}
}

/*------------------------------------------------------------------------
 *
 *      Function:    TrapezoidIntegratorSub1
 *      Description: Implements a numerical timestep integrator using
 *                   the Trapezoid integration method. Integrates only
 *					 group 1 nodes.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
void TrapezoidIntegratorSub1(Home_t *home, Node_t ***seglist1, int seglist1size)
{
        int     i, convergent, maxIterations, incrDelta;
        int     iter, globalIterError, mobIterError ,tmp = 0;;
        int     dumpErrorData = 1, doAll = 1;
		real8   errMax, globalErrMax;
        real8   oldDT, newDT;
        real8   oldx, oldy, oldz;
        real8   localVals[2], globalVals[2];
		real8	vtmp[3], currvtmp[3], vout[3], currvout[3], drloc[3], dr[3];
		real8	r[3], rold[3], r_loc[3], rold_loc[3];
		real8 	Q[3][3], Qinv[3][3];
		real8	errMaxold, xmax, ymax, zmax;
#ifdef _IMPLICIT
		real8	errx, erry, errz, thiserr, drn, relerrMax;
		Node_t	*nodemax;
#endif
        Node_t  *node;
        Param_t *param;

        param = home->param;

        oldDT = param->deltaTTsub;
        newDT = MIN(param->maxDT, param->nextDTsub);
        if (newDT <= 0.0) newDT = param->maxDT;
        param->deltaTTsub = newDT;

/*
 *      Preserve certain items of nodal data for later use.
 */
        PreserveNodalDataSub1(home, NODE_POSITION | NODE_CURR_VEL);

/*
 *		Determine the coordinate system and transformation matrix for all
 *		nodes.
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
	        node = home->nodeKeys[i];
	        if (node == (Node_t *)NULL) continue;
			node->flags |= UPDATE_NODE_J;
		}
		//AssignCoordSysSub(home, &tmp, 1);
		AssignCoordSys(home, &tmp, 1);
/*
 *      Loop until we converge on a time step.  First step is to
 *      use the current positions and velocities andreposition the
 *      nodes to where they would be after the suggested delta time.
 */
        convergent = 0;
        maxIterations = 2;
        incrDelta = 1;

        while (!convergent) {

/*
 *          Advance all nodes from their previous positions to new positions
 *          based on their current velocities and the delta T being tested.
 *          This includes all ghost nodes, so nodal velocities must previously
 *          have been distributed to neighboring domains.
 */
            AdvanceAllNodesSub1(home, oldDT);

            mobIterError = 0;
            globalVals[0] = 0.0;
            globalVals[1] = 0.0;

            for (iter = 0; iter < maxIterations; iter++) {
/*
 *              Recalculate nodal force and velocity at the new positions
 */
				FlagGroup1Nodes(home);
				//NodeForce(home, PARTIAL);
				NodeForceList(home, seglist1, seglist1size, PARTIAL );
                mobIterError = CalcNodeVelocities(home, 0, !doAll);
                CommSendVelocity(home);

/*
 *              If the mobility function was unable to iterate to
 *              converge on a velocity for one or more nodes, we
 *              just want to cut the timestep, so set the starting
 *              loop index below so we don't even bother calculating
 *              positioning errors.
 */
                if (mobIterError != 0) {
                    i = home->newNodeKeyPtr;
                } else {
                    i = 0;
                }

                errMax = 0.0;
#ifdef _IMPLICIT
				relerrMax = 0.0;
#endif

                for (/* i initializex above */; i < home->newNodeKeyPtr; i++) {

                    node = home->nodeKeys[i];
                    if (node == (Node_t *)NULL) continue;

/*
 *					If this node does not belong to group 1, skip it!.
 */
					if (node->subgroup!=1) {
						continue;
					}
                    oldx = node->oldx;
                    oldy = node->oldy;
                    oldz = node->oldz;

                    PBCPOSITION(param, node->x, node->y, node->z,
                                &oldx, &oldy, &oldz);
			
					errMaxold = errMax;

					UnpackQMatrices(node, Q, Qinv);

					r[0] = node->x; 	
					r[1] = node->y; 	
					r[2] = node->z;
					rold[0] = oldx; 	
					rold[1] = oldy; 	
					rold[2] = oldz;
					vtmp[0] = node->vX;	
					vtmp[1] = node->vY;	
					vtmp[2] = node->vZ;					
					currvtmp[0] = node->currvX;
					currvtmp[1] = node->currvY;
					currvtmp[2] = node->currvZ;

					Matrix33Vector3Multiply(Qinv, r, r_loc);
					Matrix33Vector3Multiply(Qinv, rold, rold_loc);
					Matrix33Vector3Multiply(Qinv, vtmp, vout);
					Matrix33Vector3Multiply(Qinv, currvtmp, currvout);

					if (node->ndof>0) {
						errx = r_loc[0]-rold_loc[0]-((vout[0]+currvout[0])*0.5*newDT);
					} else {
						errx = 0;
					}
					if (node->ndof>1) {
						erry = r_loc[1]-rold_loc[1]-((vout[1]+currvout[1])*0.5*newDT);
					} else {
						erry = 0;
					}
					if (node->ndof>2) {
						errz = r_loc[2]-rold_loc[2]-((vout[2]+currvout[2])*0.5*newDT);
					} else {
						errz = 0;
					}					
/*
#ifndef _IMPLICIT
                    errMax = MAX(errMax, fabs(node->x - oldx -
                                         ((node->vX+node->currvX)*0.5*newDT)));
                    errMax = MAX(errMax, fabs(node->y - oldy -
                                         ((node->vY+node->currvY)*0.5*newDT)));
                    errMax = MAX(errMax, fabs(node->z - oldz -
                                         ((node->vZ+node->currvZ)*0.5*newDT)));
#else

					errx = node->x - oldx - ((node->vX+node->currvX)*0.5*newDT);
					erry = node->y - oldy - ((node->vY+node->currvY)*0.5*newDT);
					errz = node->z - oldz - ((node->vZ+node->currvZ)*0.5*newDT);
*/
					thiserr = sqrt(errx*errx+erry*erry+errz*errz);
					drn = sqrt(pow(node->x-oldx,2)+pow(node->y-oldy,2)+pow(node->z-oldz,2));
					//xx = node->x - oldx;
					//yy = node->y - oldy;
					//zz = node->z - oldz;
					//drn = sqrt( xx*xx + yy*yy + zz*zz );
					
					errMax = MAX(errMax, thiserr);
					real8 relerrMaxold = relerrMax;
					if (thiserr>param->rTolth) {
						if (drn>param->rTolth/param->rTolrel) {
							relerrMax = MAX(relerrMax,thiserr/drn);
						} else {
							relerrMax = 2*param->rTolrel;
						}
					}
					if (relerrMax!=relerrMaxold) {
						nodemax = node;
					}

//#endif
					//if (errMax!=errMaxold) {
					//	xmax = node->oldx;
					//	ymax = node->oldy;
					//	zmax = node->oldz;
					//}				
                }
	
/*
 *              Need to find largest errMax from among all domains.
 */
#if PARALLEL
                localVals[0] = errMax;
                localVals[1] = (real8)mobIterError;

                MPI_Allreduce(localVals, globalVals, 2, MPI_DOUBLE, MPI_MAX, 
                              MPI_COMM_WORLD);

                globalErrMax = globalVals[0];;
                globalIterError = globalVals[1];
#else
                globalErrMax = errMax;
                globalIterError = mobIterError;
#endif

				//printf("subiter=%i with error=%e at %e %e %e\n",iter,
								//globalErrMax,xmax,ymax,zmax);
				printf("subiter=%i with error=%e and dt=%e\n",iter,globalErrMax,newDT);
/*
				if (relerrMax>0) {
					printf("  with %e pos(%e %e %e) oldpos(%e %e %e)\n  v(%e %e %e) oldv(%e %e %e) F(%e, %e, %e)\n",
					relerrMax,nodemax->x,nodemax->y,nodemax->z,nodemax->oldx,nodemax->oldy,nodemax->oldz,
					nodemax->vX,nodemax->vY,nodemax->vZ,nodemax->currvX,nodemax->currvY,nodemax->currvZ,
					nodemax->fX,nodemax->fY,nodemax->fZ);
					printf("  %i %e %e %e %e %e %e %e %e %e\n",nodemax->ndof,nodemax->Q11,nodemax->Q12,nodemax->Q13,
										nodemax->Q21,nodemax->Q22,nodemax->Q23,
										nodemax->Q31,nodemax->Q32,nodemax->Q33);
					for (i=0; i<nodemax->numNbrs; i++) {
						printf("  %e %e %e",nodemax->nx[i],nodemax->ny[i],nodemax->nz[i]);
					}
					printf("\n");
				}
*/
/*
 *              If any domain encountered an error iterating inside
 *              the mobility function, just go right to cutting the
 *              timestep and trying again.
 */
                if (globalIterError) {
                    iter = maxIterations;
                    continue;
                }

/*
 *              If the error is within the tolerance, we've reached
 *              convergence so we can accept this deltaT.  Otherwise
 *              reposition the nodes and try again.  Note: we need to
 *              reposition both local nodes and ghost nodes!
 */
#ifndef _IMPLICIT
                if (globalErrMax < param->rTol) {
#else
				if (globalErrMax < param->rTol && relerrMax<param->rTolrel) {
#endif
                    convergent = 1;
                    break;
                } else {
                    incrDelta = 0;
                    if (iter == maxIterations-1) {
                        continue;
                    }
                    for (i = 0; i < home->newNodeKeyPtr; i++) {
                        node = home->nodeKeys[i];
                        if (node == (Node_t *)NULL) continue;
/*
 *						If this node does not belong to group 1, don't touch it!.
 */
						if (node->subgroup!=1) {
							continue;
						}
                        oldx = node->oldx;
                        oldy = node->oldy;
                        oldz = node->oldz;

                        PBCPOSITION(param, node->x, node->y, node->z,
                                    &oldx, &oldy, &oldz);

/*
 *						Transform the velocities into 
 *						local coordinates.
 */
		                UnpackQMatrices(node, Q, Qinv);

						vtmp[0] = node->vX;
						vtmp[1] = node->vY;
						vtmp[2] = node->vZ;					
						currvtmp[0] = node->currvX;
						currvtmp[1] = node->currvY;
						currvtmp[2] = node->currvZ;

						Matrix33Vector3Multiply(Qinv, vtmp, vout);
						Matrix33Vector3Multiply(Qinv, currvtmp, currvout);

/*
 *						Calculate the updates as appropriate for the local
 *						coordinates that correspond to degrees of freedom.
 */
						if (node->ndof>0) {
							drloc[0] = (vout[0]+currvout[0])*0.5*newDT;
						} else {
							drloc[0] = 0;
						}
						if (node->ndof>1) {
							drloc[1] = (vout[1]+currvout[1])*0.5*newDT;
						} else {
							drloc[1] = 0;
						}
						if (node->ndof>2) {
							drloc[2] = (vout[2]+currvout[2])*0.5*newDT;
						} else {
							drloc[2] = 0;
						}

/*						
 *						Transform the updates into global coordinates and 
 *						modify the position.
 */
						Matrix33Vector3Multiply(Q, drloc, dr);
						
                        node->x -= node->x - oldx - dr[0];
                        node->y -= node->y - oldy - dr[1];
                        node->z -= node->z - oldz - dr[2];

                        FoldBox(param, &node->x, &node->y, &node->z);
                    }

                    node = home->ghostNodeQ;

                    while (node) {

                        oldx = node->oldx;
                        oldy = node->oldy;
                        oldz = node->oldz;

                        PBCPOSITION(param, node->x, node->y, node->z,
                                    &oldx, &oldy, &oldz);

                        node->x -= node->x - oldx -
                                   ((node->vX+node->currvX)*0.5*newDT);
                        node->y -= node->y - oldy -
                                   ((node->vY+node->currvY)*0.5*newDT);
                        node->z -= node->z - oldz -
                                   ((node->vZ+node->currvZ)*0.5*newDT);

                        FoldBox(param, &node->x, &node->y, &node->z);

                        node = node->next;
                    }

                }  /* not convergent */
            }  /* for (iter = 0; ...) */

/*
 *          If there is convergence, we've got a good delta T, otherwise
 *          cut the delta T by a configured factor and try again.
 */
            if (!convergent) {

		//if (param->deltaTTsub<=param->realdt/1000) {
		//	convergent = 1;
		//	break;
		//}

                newDT *= param->dtDecrementFact;
                param->deltaTTsub = newDT;

                if ((newDT < 1.0e-20) && (home->myDomain == 0)) {
                    Fatal("TrapezoidIntegratorSub1(): Timestep has dropped below\n"
                          "minimal threshold to %e.  Aborting!", newDT);
                }
            }

        }  /* while (!convergent) */

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
        param->deltaTTsub   = newDT;
        param->realdtsub    = newDT;

        if (incrDelta) {
            if (param->dtVariableAdjustment) {
                real8 tmp1, tmp2, tmp3, tmp4, factor;
                tmp1 = pow(param->dtIncrementFact, param->dtExponent);
                tmp2 = globalErrMax/param->rTol;
                tmp3 = 1.0 / param->dtExponent;
                tmp4 = pow(1.0/(1.0+(tmp1-1.0)*tmp2), tmp3);
                factor = param->dtIncrementFact * tmp4;
                param->nextDTsub = MIN(param->maxDT, newDT*factor);
            } else {
                param->nextDTsub = MIN(param->maxDT, newDT*param->dtIncrementFact);
            }
        } else {
            param->nextDTsub = newDT;
        }

/*
 *      Copy the nodal velocities that existed on entry to the timestep
 *      integrator.
 */
        PreserveNodalDataSub1(home, NODE_OLD_VEL);

#ifdef _FEM
/*
 *      If we're using the FEM code for simulating free surfaces, we
 *      have to handle any nodes/segments that have moved onto or past
 *      the simulation surfaces.
 *
 *      FIX ME!  Do we need to again recalculate forces/velocities for
 *               the nodes that get modified in AdjustNodePosition(), or
 *               can we live with the results until the next cycle???
 */
        AdjustNodePosition(home, 1);
#endif
        return;
}
