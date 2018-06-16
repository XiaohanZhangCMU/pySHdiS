/**************************************************************************
 *
 *      Module:      RKFIntegrator.c
 *      Description: Implements a numerical timestep integrator using
 *                   the Runge-Kutta-Fehlberg integration method, a 4th
 *                   order explicit method.
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
static void PreserveNodalData(Home_t *home, int reqType)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

			if ( reqType >= GROUP0 && node->subgroup == 0) continue;
 
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

			if ( reqType >= GROUP0 && node->subgroup == 0) {    
                node = node->next;
                continue;
            }

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
 *      Function:    RKFStep
 *      Description: Function for taking a single step of the RKF method.
 *                   This function repositions the nodes but does not
 *                   update the forces and velocities.
 *
 *-----------------------------------------------------------------------*/

void RKFStep(Home_t *home, int step, int reqType) {

        int     i;
        real8   currDT, x, y, z, f1, f2, f3, f4, f5, f6; 
        Node_t  *node;

/*
 *      Array of coefficients used by the method.
 */
        static real8 f[6][6] = {
            {1.0/4      ,  0.0        ,  0.0         ,  0.0          ,  0.0    , 0.0   },
            {3.0/32     ,  9.0/32     ,  0.0         ,  0.0          ,  0.0    , 0.0   },
            {1932.0/2197, -7200.0/2197,  7296.0/2197 ,  0.0          ,  0.0    , 0.0   },
            {439.0/216  , -8.0        ,  3680.0/513  , -845.0/4104   ,  0.0    , 0.0   },
            {-8.0/27    ,  2.0        , -3544.0/2565 ,  1859.0/4104  , -11.0/40, 0.0   },
            {16.0/135   ,  0.0        ,  6656.0/12825,  28561.0/56430, -9.0/50 , 2.0/55} 
        };
		
		f1 = f[step][0];
		f2 = f[step][1];
		f3 = f[step][2];
		f4 = f[step][3];
		f5 = f[step][4];
		f6 = f[step][5];
		
		if (reqType == FULL || reqType == GROUP0) currDT = home->param->deltaTT;
		else if              ( reqType == GROUP1) currDT = home->param->deltaTTsub;
		else if              ( reqType == GROUP2) currDT = home->param->deltaTTsub2;
		else if              ( reqType == GROUP3) currDT = home->param->deltaTTsub3;
		else if              ( reqType == GROUP4) currDT = home->param->deltaTTsub4;

        for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

			if ( reqType >= GROUP0 && node->subgroup == 0) continue;

            if (step<5) {
				node->RKFx[step] = node->vX;
				node->RKFy[step] = node->vY;
				node->RKFz[step] = node->vZ;
            }
			
			x = node->oldx + currDT*(f1*node->RKFx[0] + f2*node->RKFx[1] + f3*node->RKFx[2] + 
			                         f4*node->RKFx[3] + f5*node->RKFx[4] + f6*node->RKFx[5]);
			y = node->oldy + currDT*(f1*node->RKFy[0] + f2*node->RKFy[1] + f3*node->RKFy[2] + 
			                         f4*node->RKFy[3] + f5*node->RKFy[4] + f6*node->RKFy[5]);
			z = node->oldz + currDT*(f1*node->RKFz[0] + f2*node->RKFz[1] + f3*node->RKFz[2] + 
			                         f4*node->RKFz[3] + f5*node->RKFz[4] + f6*node->RKFz[5]);

			FoldBox(home->param,&x,&y,&z);

			node->x = x;
			node->y = y;
			node->z = z;
		}

        node = home->ghostNodeQ;

        while (node != (Node_t *)NULL) {

			if ( reqType >= GROUP0 && node->subgroup == 0) {    
                node = node->next;
                continue;
            }

			if (step<5) {
				node->RKFx[step] = node->vX;
				node->RKFy[step] = node->vY;
				node->RKFz[step] = node->vZ;
            }
			
			x = node->oldx + currDT*(f1*node->RKFx[0] + f2*node->RKFx[1] + f3*node->RKFx[2] + 
			                         f4*node->RKFx[3] + f5*node->RKFx[4] + f6*node->RKFx[5]);
			y = node->oldy + currDT*(f1*node->RKFy[0] + f2*node->RKFy[1] + f3*node->RKFy[2] + 
			                         f4*node->RKFy[3] + f5*node->RKFy[4] + f6*node->RKFy[5]);
			z = node->oldz + currDT*(f1*node->RKFz[0] + f2*node->RKFz[1] + f3*node->RKFz[2] + 
			                         f4*node->RKFz[3] + f5*node->RKFz[4] + f6*node->RKFz[5]);

			FoldBox(home->param,&x,&y,&z);

			node->x = x;
			node->y = y;
			node->z = z;

            node = node->next;
        }

}

/*------------------------------------------------------------------------
 *
 *      Function:    RKFIntegrator
 *      Description: Implements a numerical timestep integrator using
 *                   the Runge-Kutta-Fehlberg integration method.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
void RKFIntegrator(Home_t *home, int reqType)
{
        int      i, j, convergent, incrDelta, iTry;
        int      globalIterError, mobIterError;
        int      dumpErrorData = 1, doAll = 1;
		real8    errMax, globalErrMax, globalRelErrMax;
        real8    currDT, newDT;
        real8    oldx, oldy, oldz, x, y, z, dx, dy, dz;
        real8    drn, relerrMax, rg9s, tmp_rg9s;
        real8 	 er1, er2, er3, er4, er5, er6;		
		real8 	 errx, erry, errz, errnet;
        real8    localVals[3], globalVals[3];
        Node_t   *node, *node1, *node2, *node3, *node4;
        Param_t  *param;
		Subcyc_t *subcyc;
		
		param   = home->param;
		subcyc  = home->subcyc;

        // Coefficients for error calculation
		er1 =  1.0/360   ;   er4 = -2197.0/75240;
		er2 =  0.0       ;   er5 =  1.0/50      ;
		er3 = -128.0/4275;   er6 =  2.0/55      ;
		
/*
 *		Calculate the current velocities for the given group defined by reqType
 */
//		NodeForce(home, reqType);
//		mobIterError = CalcNodeVelocities(home, 0, doAll);
//       if (reqType==FULL) {
//		    CommSendVelocity(home);
//        } else {
//            CommSendVelocitySub(home, subGroup);
//        }

/*
 *      Grab the time step.
 */
		if (reqType == FULL || reqType == GROUP0) newDT = MIN(param->maxDT, param->nextDT    );
		else if              ( reqType == GROUP1) newDT = MIN(param->maxDT, param->nextDTsub );
		else if              ( reqType == GROUP2) newDT = MIN(param->maxDT, param->nextDTsub2);
		else if              ( reqType == GROUP3) newDT = MIN(param->maxDT, param->nextDTsub3);
		else if              ( reqType == GROUP4) newDT = MIN(param->maxDT, param->nextDTsub4);
		
		if (newDT <= 0.0) {
			if (reqType == FULL || reqType == GROUP0) newDT = param->maxDT;
			else newDT = param->realdt;
		}
		
		if (reqType == FULL || reqType == GROUP0) param->deltaTT     = newDT;
		else if              ( reqType == GROUP1) param->deltaTTsub  = newDT;
		else if              ( reqType == GROUP2) param->deltaTTsub2 = newDT;
		else if              ( reqType == GROUP3) param->deltaTTsub3 = newDT;
		else if              ( reqType == GROUP4) param->deltaTTsub4 = newDT;
		
/*
 *		Initialize variables needed for automatic adjustment of group0 interactions.
 */
		if (reqType == FULL || reqType == GROUP0) {
			rg9s = MAX(MAX(MAX(param->rg1,param->rg2),param->rg3),param->rg4) * 2;
			rg9s = rg9s * rg9s;
			rg9s = 5000.0 *5000.0 ;
//Amin - what should this be set to, not a constant right????
		}
		
		if (reqType == GROUP0) {
			for (i=0; i < home->newNodeKeyPtr; i++) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				node->G0_to_G4 = 0 ;
			}
		}

/*
 *      Initialize velocity arrays to zeros.
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if ( reqType >= GROUP0 && node->subgroup == 0) continue;
			for (j=0; j<6; j++) {
				node->RKFx[j] = 0.0;
				node->RKFy[j] = 0.0;
				node->RKFz[j] = 0.0;
			}
		}

        node = home->ghostNodeQ;
        while (node != (Node_t *)NULL) {
			if ( reqType >= GROUP0 && node->subgroup == 0) {
                node = node->next;
                continue;
            }
			for (j=0; j<6; j++) {
				node->RKFx[j] = 0.0;
				node->RKFy[j] = 0.0;
				node->RKFz[j] = 0.0;
			}
            node = node->next;
        }

/*
 *      Preserve certain items of nodal data for later use.
 */
        PreserveNodalData(home, reqType);

/*
 *      Loop until we converge on a time step.  First step is to
 *      use the current positions and velocities and reposition the
 *      nodes to where they would be after the suggested delta time.
 */
        convergent =  0;
		incrDelta  =  1;
		iTry       = -1;

        while (!convergent) {
			iTry++;
			mobIterError = 0;
			globalVals[0] = 0.0;
			globalVals[1] = 0.0;
			globalVals[2] = 0.0;
			errMax = 0.0;
            relerrMax = 0.0;

/*
 *			Apply the Runge-Kutta-Fehlberg integrator one step at a time, 
 *			checking for mobility errors after each velocity calculation.
 */

            for (i=0 ; i<5 ; i++) {
                RKFStep(home, i, reqType);
				
				NodeForce(home, reqType);
				mobIterError = CalcNodeVelocities(home, 0, doAll);
                if (reqType==FULL) {
		            CommSendVelocity(home);
                } else {
                    CommSendVelocitySub(home, reqType);
                }
                if (mobIterError != 0) break;
            }

            if (mobIterError != 0) {
                i = home->newNodeKeyPtr;
            } else {
                i = 0;
            }
            
            //Calculate the error
			for (/* initialized above */; i < home->newNodeKeyPtr; i++) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

				if ( reqType >= GROUP0 && node->subgroup == 0) continue;
    
				node->RKFx[5] = node->vX;
				node->RKFy[5] = node->vY;
				node->RKFz[5] = node->vZ;               

				errx = newDT*(er1*node->RKFx[0] + er2*node->RKFx[1] + er3*node->RKFx[2] + 
				              er4*node->RKFx[3] + er5*node->RKFx[4] + er6*node->RKFx[5]);
				erry = newDT*(er1*node->RKFy[0] + er2*node->RKFy[1] + er3*node->RKFy[2] + 
				              er4*node->RKFy[3] + er5*node->RKFy[4] + er6*node->RKFy[5]);
				errz = newDT*(er1*node->RKFz[0] + er2*node->RKFz[1] + er3*node->RKFz[2] + 
				              er4*node->RKFz[3] + er5*node->RKFz[4] + er6*node->RKFz[5]);

			    errnet = sqrt(errx*errx + erry*erry + errz*errz);
			    errMax = MAX (errMax, errnet);

                oldx = node->oldx;
                oldy = node->oldy;
                oldz = node->oldz;

                PBCPOSITION(param, node->x, node->y, node->z,
                                   &oldx  , &oldy  , &oldz);

                dx = node->x - oldx;
                dy = node->y - oldy;
                dz = node->z - oldz;
                drn = sqrt(dx*dx + dy*dy + dz*dz);
                if (errnet>param->rTolth) {
                    if (drn>param->rTolth/param->rTolrel) {
                        relerrMax = MAX(relerrMax,errnet/drn);
                    } else {
                        relerrMax = 2*param->rTolrel;
                    }
              	}
				
				if (reqType == GROUP0 && iTry < param->nTry) {
					if ((errnet < param->rTol ) && 
					    (errnet < param->rTolth || errnet/drn < param->rTolrel)) {
						node->G0_to_G4 = 0;
					} else {
						node->G0_to_G4 = 1;
					}
				}
			}
                
/*
 *			Need to find largest errMax from among all domains.
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
 *			If any domain encountered an error iterating inside
 *			the mobility function, just go right to cutting the
 *			timestep and trying again.
 */
			if (globalIterError) globalErrMax = 10.0*param->rTol;

/*
 *			If the error is within the tolerance, we've reached
 *			convergence so we can accept this deltaT.  Otherwise
 *			reposition the nodes and try again.  Note: we need to
 *			reposition both local nodes and ghost nodes!
 */
			if (globalErrMax<param->rTol && globalRelErrMax<param->rTolrel) {
				//Calculate final positions
				RKFStep(home, 5, reqType);
				NodeForce(home, reqType);
				mobIterError = CalcNodeVelocities(home, 0, doAll);
                if (reqType==FULL) {
		            CommSendVelocity(home);
                } else {
                    CommSendVelocitySub(home, reqType);
                }

				convergent = 1;
			} 

/*
 *			If there is convergence, we've got a good delta T, otherwise
 *			cut the delta T by a configured factor and try again.
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
/*
 *					We need to start from the old velocities. So, first, 
 *					substitute them with the old ones.  
 */
					for (i=0; i < home->newNodeKeyPtr; i++) {
						if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
						if ( reqType >= GROUP0 && node->subgroup == 0) continue;
						node->vX = node->currvX;
						node->vY = node->currvY;
						node->vZ = node->currvZ;
					}
					node = home->ghostNodeQ;
					while (node) {		
				        if ( reqType >= GROUP0 && node->subgroup == 0) {
                            node = node->next;
                            continue;
                        }	
						node->vX = node->currvX;
						node->vY = node->currvY;
						node->vZ = node->currvZ;
						node = node->next;
					}
					
					incrDelta = 0;
					newDT    *= param->dtDecrementFact;
					if (reqType == FULL || reqType == GROUP0) param->deltaTT     = newDT;
					else if              ( reqType == GROUP1) param->deltaTTsub  = newDT;
					else if              ( reqType == GROUP2) param->deltaTTsub2 = newDT;
					else if              ( reqType == GROUP3) param->deltaTTsub3 = newDT;
					else if              ( reqType == GROUP4) param->deltaTTsub4 = newDT;

					if ((newDT < 1.0e-20) && (home->myDomain == 0)) {
						Fatal("RKFIntegrator(): Timestep has dropped below\n"
							  "minimal threshold to %e.  Aborting!", newDT);
					}
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
		
		if (reqType == FULL || reqType == GROUP0) {
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
		
		if (reqType == FULL || reqType == GROUP0) param->nextDT     = newDT;
		else if              ( reqType == GROUP1) param->nextDTsub  = newDT;
		else if              ( reqType == GROUP2) param->nextDTsub2 = newDT;
		else if              ( reqType == GROUP3) param->nextDTsub3 = newDT;
		else if              ( reqType == GROUP4) param->nextDTsub4 = newDT;

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
