/**************************************************************************
 *
 *      Module:      TrapezoidIntegratorSub.c
 *      Description: Can be used with the SubcycleIntegratorForceB 
 *                   function (force based subcycling).
 *                   Implements a numerical timestep integrator using
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

void CommSendVelocitySub(Home_t *home, int subGroup);

/*
 *      NODE POSITIONING METHODS:
 *
 *      We're testing several methods for taking an initial guess at the
 *      new nodal positions in this timestep integration module.  If none
 *      has been selected in 'makefile.setup', set a default method.
 *
 */

#define DIR_TIMESTEP_ERROR "timestep_error"
#define DIR_NODAL_TIMESTEP "timestep"


#ifdef DEBUG_TIMESTEP
/*
 *      The following is just a debug function used to dump
 *      information about the nodes that are causing the 
 *      timestep to be cut.
 */
static void DumpNode(Node_t *node)
{
        int i, j;

        if (node == (Node_t *)NULL) return;

#if 1
        printf("  node(%d,%d) arms %d,  ",
               node->myTag.domainID, node->myTag.index, node->numNbrs);
        for (i = 0; i < node->numNbrs; i++) {
            printf(" (%d,%d)", node->nbrTag[i].domainID, node->nbrTag[i].index);
        }
        printf("\n");
#endif

/*
 *      Print the nodal velocity and total node force
 */
#if 1
        printf("  node(%d,%d)     position = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->x, node->y, node->z);
        printf("  node(%d,%d) old position = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->oldx, node->oldy, node->oldz);
#endif

/*
 *      Print the nodal velocity and total node force
 */
#if 1
        printf("  node(%d,%d)    v = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->vX, node->vY, node->vZ);
#endif
#if 1
        printf("  node(%d,%d)    f = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->fX, node->fY, node->fZ);
#endif


/*
 *      Print the old nodal velocity and total node force
 */
#if 0
        printf("  node(%d,%d) oldv = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->oldvX, node->oldvY, node->oldvZ);
#endif
#if 0
        printf("  node(%d,%d) oldf = (%e %e %e)\n",
               node->myTag.domainID, node->myTag.index, 
               node->oldfX, node->oldfY, node->oldfZ);
#endif


/*
 *      Print the arm specific forces
 */
#if 0
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d %d) arm[%d]-> (%d %d) f = (%e %e %e)\n",
                   node->myTag.domainID, node->myTag.index, i,       
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->armfx[i], node->armfy[i], node->armfz[i]);
        }
#endif

/*
 *      Print the burger's vector for each arm of the node
 */
#if 0
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d %d) arm[%d]-> (%d %d) b = (%f %f %f)\n",
                   node->myTag.domainID, node->myTag.index, i,       
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->burgX[i], node->burgY[i], node->burgZ[i]);
        }
#endif


/*
 *      Print the glide plane normal for each arm of the node
 */
#if 0
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d %d) arm[%d]-> (%d %d) n = (%f %f %f)\n",
                   node->myTag.domainID,node->myTag.index, i,       
                   node->nbrTag[i].domainID,node->nbrTag[i].index,
                   node->nx[i],node->ny[i],node->nz[i]);
        }
#endif

        return;
}
#endif  /* if DEBUG_TIMESTEP */


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
static void PreserveNodalData(Home_t *home, int items)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
 
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
				
				node->oldvX = node->currvX;
				node->oldvY = node->currvY;
				node->oldvZ = node->currvZ;
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
				
				node->oldvX = node->currvX;
				node->oldvY = node->currvY;
				node->oldvZ = node->currvZ;
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
 *	Function:	AdvanceAllNodes
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
static void AdvanceAllNodes(Home_t *home, int reqType)
{
		int       i;
		real8     x, y, z;
		real8     currDT;
		Node_t   *node;
		Param_t  *param;

		param = home->param;
		
		if (reqType == GROUP0) {
			param->realdt     = param->deltaTT;
			currDT            = param->realdt;
		} else if (reqType == GROUP1) {
			param->realdtsub  = param->deltaTTsub;
			currDT            = param->realdtsub ;
		} else if (reqType == GROUP2) {
			param->realdtsub2 = param->deltaTTsub2;
			currDT            = param->realdtsub2 ;
		} else if (reqType == GROUP3) {
			param->realdtsub3 = param->deltaTTsub3;
			currDT            = param->realdtsub3 ;
		} else if (reqType == GROUP4) {
			param->realdtsub4 = param->deltaTTsub4;
			currDT            = param->realdtsub4 ;
		}

		for (i = 0; i < home->newNodeKeyPtr; i++) {

			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if ( reqType >= GROUP0 && node->subgroup == 0) continue;

/*
 *			If we don't have a value for the previous velocity, assume
 *			previous velocity was same as current velocity.
 */
			if ((node->oldvX == 0.0) && (node->oldvY == 0.0) &&
				(node->oldvZ == 0.0)) {
				node->oldvX = node->currvX;
				node->oldvY = node->currvY;
				node->oldvZ = node->currvZ;
			}
			
			x = node->oldx + 0.5 * (node->currvX + node->oldvX) * currDT;
			y = node->oldy + 0.5 * (node->currvY + node->oldvY) * currDT;
			z = node->oldz + 0.5 * (node->currvZ + node->oldvZ) * currDT;

			FoldBox(param, &x, &y, &z);
			
			node->x = x;
			node->y = y;
			node->z = z;
		}

/*
 *		Also need to move Ghost nodes
 */
		node = home->ghostNodeQ;

		while (node) {

			if ( reqType >= GROUP0 && node->subgroup == 0) {    
                node = node->next;
                continue;
            }
		
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
 *      Function:    TrapezoidIntegratorSub
 *      Description: Implements a numerical timestep integrator using
 *                   the Trapezoid integration method.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
void TrapezoidIntegratorSub(Home_t *home , int reqType)
{
		int      i, convergent, maxIterations, incrDelta;
		int      iter, globalIterError, mobIterError;
		int      dumpErrorData = 1, doAll = 1;
		real8    errMax, globalErrMax, globalRelErrMax;
		real8    newDT , rg9s , tmp_rg9s;
		real8    oldx, oldy, oldz;
		real8    localVals[3], globalVals[3];
#ifdef _IMPLICIT
		real8	 errx  , erry  , errz  , thiserr  , drn;
		real8    errdX , errdY , errdZ , relerrMax;
#endif
		Node_t   *node , *node1 , *node2 , *node3 , *node4;
		Param_t  *param;
		Subcyc_t *subcyc;
#ifdef DEBUG_TIMESTEP
		real8    oldErrMax;
		Node_t   *tmpNode;
#endif

		param  = home->param;
		subcyc = home->subcyc;
		
		if (reqType == GROUP0) {
			//rg9s  = param->rg4 * param->rg4 * 4.0;
			rg9s = 5000.0 *5000.0;
			
			for (i=0; i < home->newNodeKeyPtr; i++) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				node->G0_to_G4 = 0 ;
			}
		}
		
		if (reqType == FULL || reqType == GROUP0) newDT = MIN(param->maxDT, param->nextDT    );
		else if              ( reqType == GROUP1) newDT = MIN(param->maxDT, param->nextDTsub );
		else if              ( reqType == GROUP2) newDT = MIN(param->maxDT, param->nextDTsub2);
		else if              ( reqType == GROUP3) newDT = MIN(param->maxDT, param->nextDTsub3);
		else if              ( reqType == GROUP4) newDT = MIN(param->maxDT, param->nextDTsub4);
		
		if (newDT <= 0.0) newDT = param->maxDT;
		
		if (reqType == FULL || reqType == GROUP0) param->deltaTT     = newDT;
		else if              ( reqType == GROUP1) param->deltaTTsub  = newDT;
		else if              ( reqType == GROUP2) param->deltaTTsub2 = newDT;
		else if              ( reqType == GROUP3) param->deltaTTsub3 = newDT;
		else if              ( reqType == GROUP4) param->deltaTTsub4 = newDT;
		
/*
 *      Preserve certain items of nodal data for later use.
 */
        PreserveNodalData(home, NODE_POSITION | NODE_CURR_VEL);
		
/*
 *      Loop until we converge on a time step.  First step is to
 *      use the current positions and velocities andreposition the
 *      nodes to where they would be after the suggested delta time.
 */
        convergent = 0;
        maxIterations = 2;
        incrDelta = 1;
		
		int iTry = -1;
		real8  newDT0 = newDT;
		
        while (!convergent) {
			iTry++;

/*
 *          Advance all nodes from their previous positions to new positions
 *          based on their current velocities and the delta T being tested.
 *          This includes all ghost nodes, so nodal velocities must previously
 *          have been distributed to neighboring domains.
 */
            AdvanceAllNodes(home, reqType);

            mobIterError = 0;
            globalVals[0] = 0.0;
            globalVals[1] = 0.0;
            globalVals[2] = 0.0;

            for (iter = 0; iter < maxIterations; iter++) {
/*
 *              Recalculate nodal force and velocity at the new positions
 */
				NodeForce(home, reqType);
				mobIterError = CalcNodeVelocities(home, 0, doAll);
                if (reqType==FULL) {
		            CommSendVelocity(home);
                } else {
                    CommSendVelocitySub(home, reqType);
                }

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
                    if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
					if (reqType >= GROUP0 && node->subgroup == 0) continue;

                    oldx = node->oldx;
                    oldy = node->oldy;
                    oldz = node->oldz;

                    PBCPOSITION(param, node->x, node->y, node->z,
                                &oldx, &oldy, &oldz);

#ifdef DEBUG_TIMESTEP
                    oldErrMax = errMax;
#endif

					errx = node->x - oldx - ((node->vX+node->currvX)*0.5*newDT);
					erry = node->y - oldy - ((node->vY+node->currvY)*0.5*newDT);
					errz = node->z - oldz - ((node->vZ+node->currvZ)*0.5*newDT);
					
#ifndef _IMPLICIT
                    errMax = MAX(errMax, fabs(errx));
                    errMax = MAX(errMax, fabs(erry));
                    errMax = MAX(errMax, fabs(errz));
#else
					thiserr = sqrt( errx*errx + erry*erry + errz*errz );
					errMax  = MAX ( errMax , thiserr);
					
					errdX = node->x - oldx;
					errdY = node->y - oldy;
					errdZ = node->z - oldz;
					drn   = sqrt( errdX*errdX + errdY*errdY + errdZ*errdZ );
					
					if (thiserr > param->rTolth ) {
						if (drn > param->rTolth / param->rTolrel ) {
							relerrMax = MAX( relerrMax , thiserr/drn );
						} else {
							relerrMax = 2.0*param->rTolrel;
						}
					}
					
					if (reqType == GROUP0 && iTry < param->nTry) {
						if ((thiserr < param->rTol ) && 
						    (thiserr < param->rTolth || thiserr/drn < param->rTolrel)) {
							node->G0_to_G4 = 0;
						} else {
							node->G0_to_G4 = 1;
						}
					}
#endif

#ifdef DEBUG_TIMESTEP
					if (errMax > oldErrMax) tmpNode = node;
#endif
                }
/*
 *              Need to find largest errMax from among all domains.
 */
#if PARALLEL
				localVals[0] = errMax;
				localVals[1] = relerrMax;
				localVals[2] = (real8)mobIterError;

				MPI_Allreduce(localVals, globalVals, 3, MPI_DOUBLE, MPI_MAX, 
							  MPI_COMM_WORLD);

				globalErrMax    = globalVals[0];
				globalRelErrMax = globalVals[1];
				globalIterError = globalVals[2];
#else
				globalErrMax    = errMax;
				globalRelErrMax = relerrMax;
				globalIterError = mobIterError;
#endif

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
				if (globalErrMax < param->rTol && globalRelErrMax < param->rTolrel) {
#endif
				    NodeForce(home, reqType);
				    mobIterError = CalcNodeVelocities(home, 0, doAll);
                    if (reqType==FULL) {
		                CommSendVelocity(home);
                    } else {
                        CommSendVelocitySub(home, reqType);
                    }
                    convergent = 1;
                    break;
                } else {
                    if (reqType == GROUP0) {
						if (iTry >= param->nTry) incrDelta = 0;
					} else {
						incrDelta = 0;
					}
					
                    if (iter == maxIterations-1) continue;
                    
                    for (i = 0; i < home->newNodeKeyPtr; i++) {
                        if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
						if ( reqType >= GROUP0 && node->subgroup == 0) continue;

                        oldx = node->oldx;
                        oldy = node->oldy;
                        oldz = node->oldz;

                        PBCPOSITION(param, node->x, node->y, node->z,
                                    &oldx, &oldy, &oldz);

                        node->x = oldx + ((node->vX+node->currvX)*0.5*newDT);
                        node->y = oldy + ((node->vY+node->currvY)*0.5*newDT);
                        node->z = oldz + ((node->vZ+node->currvZ)*0.5*newDT);

                        FoldBox(param, &node->x, &node->y, &node->z);
                    }

                    node = home->ghostNodeQ;

                    while (node) {

						if ( reqType >= GROUP0 && node->subgroup == 0) {
                            node = node->next;
                            continue;
                        }

                        oldx = node->oldx;
                        oldy = node->oldy;
                        oldz = node->oldz;

                        PBCPOSITION(param, node->x, node->y, node->z,
                                    &oldx, &oldy, &oldz);

                        node->x = oldx + ((node->vX+node->currvX)*0.5*newDT);
                        node->y = oldy + ((node->vY+node->currvY)*0.5*newDT);
                        node->z = oldz + ((node->vZ+node->currvZ)*0.5*newDT);

                        FoldBox(param, &node->x, &node->y, &node->z);

                        node = node->next;
                    }

                }  /* not convergent */
            }  /* for (iter = 0; ...) */

#ifdef DEBUG_TIMESTEP_ERROR
            if (dumpErrorData) {
                DumpTimestepError(home, newDT);
                dumpErrorData = 0;
            }
#endif

/*
 *          If there is convergence, we've got a good delta T, otherwise
 *          cut the delta T by a configured factor and try again.
 */
            if (!convergent) {
				if ((iTry < param->nTry) && (reqType == GROUP0)) {
					tmp_rg9s = rg9s * (iTry+1) * (iTry+1);
					
					for (i = 0; i < subcyc->SegSegListG0_cnt ; i++) {
						if (subcyc->SegSegListG0[i].flag == 0) continue;
						
						node1 = subcyc->SegSegListG0[i].seg1->node1;
						node2 = subcyc->SegSegListG0[i].seg1->node2;
						node3 = subcyc->SegSegListG0[i].seg2->node1;
						node4 = subcyc->SegSegListG0[i].seg2->node2;
						
						if ((node1->G0_to_G4 == 1) || (node2->G0_to_G4 == 1) ||
						    (node3->G0_to_G4 == 1) || (node4->G0_to_G4 == 1)) {
							
							if (subcyc->SegSegListG0[i].dist2 > tmp_rg9s) continue;
							subcyc->SegSegListG0[i].flag = 0;
						
						
							if (subcyc->SegSegListG4 == NULL || 
							    subcyc->SegSegListG4_cnt >= subcyc->SegSegListG4_siz) {
								subcyc->SegSegListG4_siz += 100;
								subcyc->SegSegListG4 = (SegSeg_t *)realloc(subcyc->SegSegListG4, 
														sizeof(SegSeg_t) * subcyc->SegSegListG4_siz);
							}
							
							subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].seg1  = subcyc->SegSegListG0[i].seg1;
							subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].seg2  = subcyc->SegSegListG0[i].seg2;
							subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].flag  = 1;
							subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].setSeg1Forces = 1;
							subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].setSeg2Forces = 1;
							subcyc->SegSegListG4_cnt++;
						}
					}

				//	FlagSubcycleNodes(home,GROUP0);
					
				} else {
					newDT *= param->dtDecrementFact;
						
					if      (reqType == GROUP0) param->deltaTT     = newDT;
					else if (reqType == GROUP1) param->deltaTTsub  = newDT;
					else if (reqType == GROUP2) param->deltaTTsub2 = newDT;
					else if (reqType == GROUP3) param->deltaTTsub3 = newDT;
					else if (reqType == GROUP4) param->deltaTTsub4 = newDT;
						
					if ((newDT < 1.0e-20) && (home->myDomain == 0)) {
						Fatal("TrapezoidIntegratorSub(): Timestep has dropped below\n"
							  "minimal threshold to %e.  Aborting!", newDT);
					}
				

#ifdef DEBUG_TIMESTEP
					if ((home->myDomain == 0) && (globalIterError)) {
						printf(" +++ Cut timestep to %e for mobility "
							   "non-convergence\n", newDT);
					}
/*
 *              	If this is this domain with the node causing timestep to drop
 *              	dump some info on the node limiting the timestep.
 */
					if ((globalErrMax == errMax) && (globalIterError == 0)) {
						printf("  Cut timestep for (%d,%d):  errMax %e, newDT %e\n",
							   tmpNode->myTag.domainID, tmpNode->myTag.index,
							   errMax, newDT);
						DumpNode(tmpNode);
					}
#endif
				}
			}

        }  /* while (!convergent) */

#ifdef DEBUG_NODAL_TIMESTEP
        DumpPerNodeTimestep(home, newDT);
#endif

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
                real8 tmp1, tmp2, tmp3, tmp4, factor;
                tmp1 = pow(param->dtIncrementFact, param->dtExponent);
                tmp2 = globalErrMax/param->rTol;
                tmp3 = 1.0 / param->dtExponent;
                tmp4 = pow(1.0/(1.0+(tmp1-1.0)*tmp2), tmp3);
                factor = param->dtIncrementFact * tmp4;
				
				newDT = MIN(param->maxDT, newDT*factor);
            } else {
				newDT = MIN(param->maxDT, newDT*param->dtIncrementFact);
            }
        }
		
		if      ( reqType == GROUP0) param->nextDT     = newDT;
		else if ( reqType == GROUP1) param->nextDTsub  = newDT;
		else if ( reqType == GROUP2) param->nextDTsub2 = newDT;
		else if ( reqType == GROUP3) param->nextDTsub3 = newDT;
		else if ( reqType == GROUP4) param->nextDTsub4 = newDT;

/*
 *      Copy the nodal velocities that existed on entry to the timestep
 *      integrator.
 */
        PreserveNodalData(home, NODE_OLD_VEL);
		
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
