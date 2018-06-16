
#include "Home.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     RemoteSigbSubForce
 * 
 *      Description:  This subroutine calculates the force from the sigbFMM
 *                    contribution pre-computed on a segment of group 1. 
 *
 *-------------------------------------------------------------------------*/
void RemoteSigbSubForce(Home_t *home, int SegID, double p1f[3], double p2f[3])
{
		int      i, numPoints;
		real8    *positions, *weights;
		real8    temp, mult1, mult2;
		real8    sigbx, sigby, sigbz;
		real8    fLinvx, fLinvy, fLinvz;
		real8    p1[3], p2[3];
		real8    pspanx, pspany, pspanz;
		Param_t  *param;
		Subcyc_t *subcyc;
		Node_t   *node1, *node2;
		
		param  = home->param;
		subcyc = home->subcyc;
		
		numPoints = param->fmNumPoints;
		positions = home->glPositions;
		weights   = home->glWeights;
		
		p1f[0] = 0.0;
		p1f[1] = 0.0;
		p1f[2] = 0.0;

		p2f[0] = 0.0;
		p2f[1] = 0.0;
		p2f[2] = 0.0;
		
		node1 = subcyc->SegListG1[SegID].seg->node1;
		node2 = subcyc->SegListG1[SegID].seg->node2;
		
		p1[X] = node1->x;
		p1[Y] = node1->y;
		p1[Z] = node1->z;

		p2[X] = node2->x;
		p2[Y] = node2->y;
		p2[Z] = node2->z;
		
		PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);
		
		pspanx = 0.5 * (p2[X]-p1[X]);
		pspany = 0.5 * (p2[Y]-p1[Y]);
		pspanz = 0.5 * (p2[Z]-p1[Z]);
	
		for (i = 0; i < numPoints; i++) {
		
			sigbx = subcyc->sigbFMM[SegID][i*3+0];
			sigby = subcyc->sigbFMM[SegID][i*3+1];
			sigbz = subcyc->sigbFMM[SegID][i*3+2];

			fLinvx = (sigby*pspanz-sigbz*pspany);
			fLinvy = (sigbz*pspanx-sigbx*pspanz);
			fLinvz = (sigbx*pspany-sigby*pspanx);

			temp = weights[i]*positions[i];
			mult1 = weights[i]+temp;

			p2f[0] = p2f[0] + fLinvx*mult1;
			p2f[1] = p2f[1] + fLinvy*mult1;
			p2f[2] = p2f[2] + fLinvz*mult1;

			mult2 = weights[i]-temp;

			p1f[0] = p1f[0] + fLinvx*mult2;
			p1f[1] = p1f[1] + fLinvy*mult2;
			p1f[2] = p1f[2] + fLinvz*mult2;    
		}
	
		return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     RemoteSigbOneSegPos
 *      Description:  This subroutine calculates the stress on a single
 *                    segment from all segments in remote cells (FMM). 
 *                    The stress is evaluated at a position specified 
 *                    along the segment.
 * 					
 *                    IMPORTANT: this function requires that the <node1>
 *                    argument corresponds to the node owning the segment!
 *
 *      Arguments:
 *          node1   Pointer to node owning the segment
 *          node2   Pointer to 2nd node of segment
 *          pos     Position in [-1,1] at wich the stress is evaluated
 *          sigb    Array in which to return the remote stress components
 *
 *-------------------------------------------------------------------------*/
void RemoteSigbOneSegPos(Home_t *home, Node_t *node1, Node_t *node2,
                         real8 pos, real8 sigb[3])
{
        int       i, j;
        int       cx, cy, cz, cellID, armID;
        real8     p1[3], p2[3], R[3], evalPos[3];
        real8     burg[3], sigma[3][3];
        real8     pmidx, pmidy, pmidz;
        real8     pspanx, pspany, pspanz;
        Param_t   *param;
        FMLayer_t *layer;
        Cell_t    *cell;
        FMCell_t  *FMcell;
        
        param = home->param;
        
        for (i = 0; i < 3; i++) sigb[i] = 0.0;
        
        if (param->fmEnabled) {
			
			layer = &home->fmLayer[param->fmNumLayers-1];
			
			armID = GetArmID(home, node1, node2);
			
			p1[X] = node1->x;
			p1[Y] = node1->y;
			p1[Z] = node1->z;

			p2[X] = node2->x;
			p2[Y] = node2->y;
			p2[Z] = node2->z;
			
			burg[X] = node1->burgX[armID];
			burg[Y] = node1->burgY[armID];
			burg[Z] = node1->burgZ[armID];

			PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);
			
/*
 *      	Find the indices (not shifted for ghost cells) of the cell
 *      	owning <node> and convert to the corresponding cellID at
 *      	the lowest FM layer.
 */
			cell = home->cellKeys[node1->cellIdx];

			cx = cell->xIndex;
			cy = cell->yIndex;
			cz = cell->zIndex;

			cx--; cy--; cz--;
			cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);
			
			FMcell = LookupFMCell(layer->cellTable, cellID);
			
/*
 *      	If PBC is enabled, the segment endpoints *may* have 
 *      	moved outside the periodic boundaries and been folded
 *      	back into the far side of the problem space.  In case
 *      	this has happened, we need to adjust the coordinates
 *      	to that of their periodic images closest to the cell
 *      	center used for the taylor expansion.
 */
			PBCPOSITION(param, FMcell->cellCtr[X], FMcell->cellCtr[Y],
						FMcell->cellCtr[Z], &p1[X], &p1[Y], &p1[Z]);
			PBCPOSITION(param, FMcell->cellCtr[X], FMcell->cellCtr[Y],
						FMcell->cellCtr[Z], &p2[X], &p2[Y], &p2[Z]);

			pmidx  = 0.5 * (p2[X]+p1[X]);
			pmidy  = 0.5 * (p2[Y]+p1[Y]);
			pmidz  = 0.5 * (p2[Z]+p1[Z]);

			pspanx = 0.5 * (p2[X]-p1[X]);
			pspany = 0.5 * (p2[Y]-p1[Y]);
			pspanz = 0.5 * (p2[Z]-p1[Z]);
			
			evalPos[X] = pmidx+pspanx*pos;
            evalPos[Y] = pmidy+pspany*pos;
            evalPos[Z] = pmidz+pspanz*pos;

			R[X] = evalPos[X] - FMcell->cellCtr[X];
			R[Y] = evalPos[Y] - FMcell->cellCtr[Y];
			R[Z] = evalPos[Z] - FMcell->cellCtr[Z];
			
			EvalTaylor(param->fmTaylorOrder, R, FMcell->taylorCoeff, sigma);
			
			
			sigb[0] += sigma[0][0]*burg[X] +
                       sigma[0][1]*burg[Y] +
                       sigma[0][2]*burg[Z];

			sigb[1] += sigma[1][0]*burg[X] +
                       sigma[1][1]*burg[Y] +
                       sigma[1][2]*burg[Z];

			sigb[2] += sigma[2][0]*burg[X] +
                       sigma[2][1]*burg[Y] +
                       sigma[2][2]*burg[Z];

		}

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     RemoteSigbSub
 * 
 *      Description:  This subroutine pre-calculates the sigb contribution
 *                    from remote cells (FMM) on each segment of group 1. 
 *                    The stress is evaluated at the Gauss points along
 *                    the segment.
 *
 *-------------------------------------------------------------------------*/
void RemoteSigbSub(Home_t *home)
{
		int      i, j, numPoints;
		Param_t  *param;
		Subcyc_t *subcyc;
		Node_t   *node1, *node2;
		real8    pos, sigb[3];

        param  = home->param;
        subcyc = home->subcyc;
        
        
        if (param->forceCutOff == 0) {
			if (param->fmEnabled) {
				
				TimerStart(home, REMOTE_FORCE);
        
/*
 * 				Zero the FMM stress contribution array
 */
				numPoints = param->fmNumPoints;
				
				for (i = 0; i < subcyc->SegListG1_cnt; i++) {
					subcyc->sigbFMM[i] = (real8*)malloc(numPoints*3*sizeof(real8));
				}

/*				Loop thru all the segments of group 1 and store 
 *              the FMM stress contribution (far-field) to be 
 *				used while subcycling group 1
 */
				for (i = 0; i < subcyc->SegListG1_cnt; i++) {
					if (subcyc->SegListG1[i].flag == 0) continue;
					
					node1 = subcyc->SegListG1[i].seg->node1;
					node2 = subcyc->SegListG1[i].seg->node2;
					
					for (j = 0; j < numPoints; j++) {
						
						pos = home->glPositions[j];
						RemoteSigbOneSegPos(home, node1, node2, pos, sigb);
						
						subcyc->sigbFMM[i][j*3+0] = sigb[0];
						subcyc->sigbFMM[i][j*3+1] = sigb[1];
						subcyc->sigbFMM[i][j*3+2] = sigb[2];
					}
					
				}

				TimerStop(home, REMOTE_FORCE);
				
			}	
		}
}

/*---------------------------------------------------------------------------
 *
 *      Function:     NodeForceList2
 *
 *-------------------------------------------------------------------------*/
void NodeForceList2(Home_t *home, SegSeg_t *SegSegList, int SegSegListSize, 
                                  Segm_t   *SegList   , int SegListSize, int subGroup)
{
		int      i, j, k;
        int      setSeg1Forces, setSeg2Forces, coreOnly;
        int      armID12, armID21, armID34, armID43;
		real8    x1, y1, z1, x2, y2, z2;
		real8    bx1, by1, bz1;
		real8    node1SegForce[3], node2SegForce[3];
        real8    a, MU, NU, Ecore;
        real8    pos1[3], pos2[3], burg[3];
        real8    dx, dy, dz ;
        real8    f1[3], f2[3], f3[3], f4[3];
        real8    sigb[3], extStress[3][3];
		real8	 eps;
        Node_t   *node1, *node2, *node3, *node4;
        Param_t  *param;
		Subcyc_t *subcyc;

        param    = home->param;
        subcyc   = home->subcyc;
        a        = param->rc;
        MU       = param->shearModulus;
        NU       = param->pois;
		coreOnly = !param->elasticinteraction;
		
/*
 *		If we are using line tension only (elasticinteraction=0), then a 
 *		different core energy is used (since the core energy is an 
 *		approximation for the total line energy rather than just that due to core 
 *		effects). Also, the tolerance for zero-length segments is different
 *		(artifact of NodeForce and LocalSegForce code). 
 */
		if (coreOnly) {
			Ecore = 0.5*param->TensionFactor*MU;
			eps = 1e-6;
		} else {			
			Ecore = param->Ecore;
			eps = 1e-20;
		}

        extStress[0][0] = param->appliedStress[0];
        extStress[1][1] = param->appliedStress[1];
        extStress[2][2] = param->appliedStress[2];
        extStress[1][2] = param->appliedStress[3];
        extStress[2][0] = param->appliedStress[4];
        extStress[0][1] = param->appliedStress[5];
        extStress[2][1] = extStress[1][2];
        extStress[0][2] = extStress[2][0];
        extStress[1][0] = extStress[0][1];


/*
 *		Reset all segment forces to zero (these forces are passed
 *		between the processors when the parallel mode is enabled).
 *		Also toggle all forcesSet to 0, they will be toggled back on 
 *		as forces are computed.
 */
		int cellNum0;
		for (cellNum0 = 0; cellNum0 < home->cellCount; cellNum0++) {
			for (i = 0; i < subcyc->totalSegCounts[cellNum0]; i++) {
				for (k = 0; k < 3; k++) {
					subcyc->cellSegLists[cellNum0][i].f1[k] = 0.0;
					subcyc->cellSegLists[cellNum0][i].f2[k] = 0.0;
					
				}
				subcyc->cellSegLists[cellNum0][i].forcesSet = 0;
			}
		}
	
/*
 *		Calculate the forces
 */
		#ifdef _OPENMP
	/*	#pragma omp parallel for default(none) schedule(dynamic,1) \*/
		#pragma omp parallel for default(none) schedule(static) \
			private(i , x1 , x2 , bx1 , dx , node1 , node3 , f1 , node1SegForce , armID12) \
			private(j , y1 , y2 , by1 , dy , node2 , node4 , f2 , node2SegForce , armID21) \
			private(k,  z1 , z2 , bz1 , dz , sigb  , f3 , f4 , armID34 ,  armID43)\
			shared (a , MU , home  , SegList    , SegListSize    , extStress) \
			shared (    NU , param , SegSegList , SegSegListSize , Ecore    )
		#endif
		for (i = 0; i < SegListSize + SegSegListSize ; i++) {
			if ( i < SegListSize ) {
				if (SegList[i].flag == 0) continue;
/*
 *				Zero out some arrays in which we'll accumulate nodal forces
 *				for the segment until we're done with the segment.  Since the
 *				AddtoArmForces() function locks the nodal info for update,
 *				we avoid calling that function until the end of the loop
 */
				VECTOR_ZERO(node1SegForce);
				VECTOR_ZERO(node2SegForce);

/*
 *				Before we calculate seg/seg forces, calculate all
 *				the forces affecting this single native segment
 *				(i.e. self force, pk force, etc) so those forces
 *				will be included when the seg/seg forces are
 *				communicated to remote domains.
 *
 *				Note: if fmm is enabled, calculate the remote force on
 *				each native segment here otherwise, the remote sigb
 *				has previously been computed and force will be calculated
 *				using that data.
 *
 *				This assumes node1 owns the segment!
 */
				node1   = SegList[i].seg->node1;
				node2   = SegList[i].seg->node2;
				armID12 = SegList[i].seg->armID12;
				armID21 = SegList[i].seg->armID21;

				x1 = node1->x;
				y1 = node1->y;
				z1 = node1->z;

				bx1 = node1->burgX[armID12];
				by1 = node1->burgY[armID12];
				bz1 = node1->burgZ[armID12];

				dx = node2->x - x1;
				dy = node2->y - y1;
				dz = node2->z - z1;

				ZImage(param, &dx, &dy, &dz);

				x2 = x1 + dx;
				y2 = y1 + dy;
				z2 = z1 + dz;
/*
 *          	Add in force due to self stress
 */
				SelfForce(0, MU, NU, bx1, by1, bz1, x1, y1, z1, x2, y2, z2,
						  a, Ecore, f1, f2);

				VECTOR_ADD(node1SegForce, f1);
				VECTOR_ADD(node2SegForce, f2);

/*
 *       	   Add in PK force from external stress
 */
				ExtPKForce(extStress, bx1, by1, bz1, x1, y1, z1,
						   x2, y2, z2, f1, f2);

				VECTOR_ADD(node1SegForce, f1);
				VECTOR_ADD(node2SegForce, f2);

/*
 *      	    If we're including osmotic forces, add those in now
 */
				if (param->vacancyConcEquilibrium > 0.0) {

					OsmoticForce(home, x1, y1, z1, x2, y2, z2,
								 bx1, by1, bz1, f1, f2);

					VECTOR_ADD(node1SegForce, f1);
					VECTOR_ADD(node2SegForce, f2);
				}

/*
 *              Add in remote force from fast-multipole method.  Arm
 *              specific forces don't need to be explicitly updated
 *              here because they are updated within RemoteForceOneSeg().
 */
                if (param->forceCutOff == 0) {
				    if (param->fmEnabled) {
						
#ifdef _PRECOMPUTE_SUB_FMM
/*
 *                      Use the pre-computed sigb FMM contribution values
 *                      to calculate the FMM force.
 * 
 */
						RemoteSigbSubForce(home, i, f1, f2);
						
						VECTOR_ADD(node1SegForce, f1);
						VECTOR_ADD(node2SegForce, f2);
#else
/*
 *                      Always pass the node owned by this domain as the 
 *                      first node for remote force calcs to make sure we 
 *                      have the taylor expansion coefficients of the FM cell.
 */					
						if (node1->myTag.domainID == home->myDomain) {
							RemoteForceOneSeg(home, node1, node2, f1, f2);
						} else {
							RemoteForceOneSeg(home, node2, node1, f2, f1);
						}
						
						for (k = 0; k < 3; k++) {
							SegList[i].seg->f1[k] += f1[k];
							SegList[i].seg->f2[k] += f2[k];
						}
#endif
				    }
                }

#if !defined _FEM & !defined _SHDIS
				if (param->fmEnabled == 0 && param->forceCutOff == 0)
#endif
				{
					sigb[0] = node1->sigbRem[armID12*3];
					sigb[1] = node1->sigbRem[armID12*3+1];
					sigb[2] = node1->sigbRem[armID12*3+2];
/*
 *					Add in PK force from non-fmm remote stress and
 *					FEM stress (if applicable).
 */
					PKForce(sigb, x1, y1, z1, x2, y2, z2, f1, f2);

          //printf("sigb[0]= %g\n", sigb[0]);
          //printf("sigb[1]= %g\n", sigb[1]);
          //printf("sigb[2]= %g\n", sigb[2]);
          //printf("f1[0]= %g\n", f1[0]);
          //printf("f1[1]= %g\n", f1[1]);
          //printf("f1[2]= %g\n", f1[2]);
          //printf("f2[0]= %g\n", f2[0]);
          //printf("f2[1]= %g\n", f2[1]);
          //printf("f2[2]= %g\n", f2[2]);


					VECTOR_ADD(node1SegForce, f1);
					VECTOR_ADD(node2SegForce, f2);
				}

/*
 *				We've accumulated various forces on the segment, so now
 *				increment the segment forces (locking will be done within
 *				the update function).
 */		
				
				#ifdef _OPENMP
					#pragma omp atomic
					node1->armfx[armID12] += node1SegForce[0];
					#pragma omp atomic
					node1->armfy[armID12] += node1SegForce[1];
					#pragma omp atomic
					node1->armfz[armID12] += node1SegForce[2];
		
					#pragma omp atomic
					node2->armfx[armID21] += node2SegForce[0];
					#pragma omp atomic
					node2->armfy[armID21] += node2SegForce[1];
					#pragma omp atomic
					node2->armfz[armID21] += node2SegForce[2];
				#else
					AddtoArmForce(node1, armID12, node1SegForce);
					AddtoArmForce(node2, armID21, node2SegForce);
				#endif

/*
 *				If we are communicating seg forces in parallel,
 *				toggle the forcesSet flag.
 */
				if (param->sendSubGroupForc)
					SegList[i].seg->forcesSet = 1;
					
				
				for (k = 0; k < 3; k++) {
					SegList[i].seg->f1[k] += node1SegForce[k];
					SegList[i].seg->f2[k] += node2SegForce[k];
				}

			}  /* end if ( i < SegListSize ) */
			else {
/*
 *				compute the force between segment pairs
 */
				j = i - SegListSize;
				if (SegSegList[j].flag == 0) continue;
				
				node1 = SegSegList[j].seg1->node1;
				node2 = SegSegList[j].seg1->node2;
				node3 = SegSegList[j].seg2->node1;
				node4 = SegSegList[j].seg2->node2;

				ComputeForces(home, node1, node2, node3, node4,
						  f1, f2, f3, f4);

				if (SegSegList[j].setSeg1Forces) {
					SegSegList[j].seg1->forcesSet = 1;

					armID12 = SegSegList[j].seg1->armID12;
					armID21 = SegSegList[j].seg1->armID21;
					
					#ifdef _OPENMP
						#pragma omp atomic
						node1->armfx[armID12] += f1[0];
						#pragma omp atomic
						node1->armfy[armID12] += f1[1];
						#pragma omp atomic
						node1->armfz[armID12] += f1[2];
			
						#pragma omp atomic
						node2->armfx[armID21] += f2[0];
						#pragma omp atomic
						node2->armfy[armID21] += f2[1];
						#pragma omp atomic
						node2->armfz[armID21] += f2[2];
					#else
						AddtoArmForce(node1, armID12, f1);
						AddtoArmForce(node2, armID21, f2);
					#endif
					
					for (k = 0; k < 3; k++) {
						SegSegList[j].seg1->f1[k] += f1[k];
						SegSegList[j].seg1->f2[k] += f2[k];
					}
				}

				if (SegSegList[j].setSeg2Forces) {
					SegSegList[j].seg2->forcesSet = 1;

					armID34 = SegSegList[j].seg2->armID12;
					armID43 = SegSegList[j].seg2->armID21;
					
					#ifdef _OPENMP
						#pragma omp atomic
						node3->armfx[armID34] += f3[0];
						#pragma omp atomic
						node3->armfy[armID34] += f3[1];
						#pragma omp atomic
						node3->armfz[armID34] += f3[2];
				
						#pragma omp atomic
						node4->armfx[armID43] += f4[0];
						#pragma omp atomic
						node4->armfy[armID43] += f4[1];
						#pragma omp atomic
						node4->armfz[armID43] += f4[2];      
					#else
						AddtoArmForce(node3, armID34, f3);
						AddtoArmForce(node4, armID43, f4);
					#endif
					
					for (k = 0; k < 3; k++) {
						SegSegList[j].seg2->f1[k] += f3[k];
						SegSegList[j].seg2->f2[k] += f4[k];
					}
				}
			}
		}
		
#ifdef PARALLEL
/*
 *      Bump up the count of segments that will be sent to
 *      any remote domain owning one of the nodes in any
 *      segment whose forces were updated by this domain.
 */
		if (subGroup <= GROUP0 || param->sendSubGroupForc) {
			int  cellNum , cellSegCnt , numDomains , homeDomain , homeCells , sendDomCnt=0;
			int        *sendDomList , *localMsgCnts , *globalMsgCnts;
			Segment_t  *segList0;

			numDomains    = home->numDomains;
			homeDomain    = home->myDomain;
			homeCells     = home->cellCount;
			sendDomList   = (int *)NULL;
			localMsgCnts  = (int *)calloc(1, sizeof(int) * numDomains);
			globalMsgCnts = (int *)calloc(1, sizeof(int) * numDomains);

			for (cellNum = 0; cellNum < homeCells; cellNum++) {
				cellSegCnt = subcyc->totalSegCounts[cellNum];
				segList0   = subcyc->cellSegLists  [cellNum];
				for (i = 0; i < cellSegCnt; i++) {
					if (segList0[i].forcesSet == 1) {
						IncrDomSegCommCnts(home, segList0[i].node1, segList0[i].node2,
										   &sendDomCnt, &sendDomList, localMsgCnts);

					}
				}
			}
		
/*
 *			Now we need to communicate the newly computed segment forces
 *			to all the appropriate remote domains.  Must NOT include this
 *			communication time with force calc time, but we do want to time
 *			the communication phase.
 */
			TimerStop(home, LOCAL_FORCE);
			TimerStop(home, CALC_FORCE);
			TimerStart(home, SEGFORCE_COMM);

			MPI_Allreduce(localMsgCnts, globalMsgCnts, numDomains,
						  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

			CommSendSegments(home, globalMsgCnts[homeDomain], sendDomCnt,
							 sendDomList, subcyc->cellSegLists, subcyc->totalSegCounts);

			TimerStop(home, SEGFORCE_COMM);
			TimerStart(home, LOCAL_FORCE);
			TimerStart(home, CALC_FORCE);
			
			if (sendDomList != (int *)NULL) free(sendDomList);
			free(globalMsgCnts);
			free(localMsgCnts);
		}
#endif
		
/*
 *      We should now have updated forces for nodes/segments
 *      so now do a quick loop through all local nodes and set
 *      the nodes' total forces to the sum of their arms' forces.
 */
		TimerStart(home, ADD_FORCE);
		
		#ifdef _OPENMP
	/*	#pragma omp parallel for default(none) schedule(dynamic,1) \ */
		#pragma omp parallel for default(none) schedule(static) \
			private(i , j , node1) \
			shared (home  )
		#endif
		for (i = 0; i < home->newNodeKeyPtr; i++) {

			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if (node1->subgroup == 0) continue;

			node1->fX = 0.0;
			node1->fY = 0.0;
			node1->fZ = 0.0;

			for (j = 0; j < node1->numNbrs; j++) {
				node1->fX += node1->armfx[j];
				node1->fY += node1->armfy[j];
				node1->fZ += node1->armfz[j];
			}
		}
		TimerStop(home, ADD_FORCE);

        return;
}
