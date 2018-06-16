
#include "Home.h"
void   SortNodes(Home_t *home, real8 *maxsep);
void ZeroLongRangeNodeForces(Home_t *home, int reqType);

void NodeForceList(Home_t *home, Node_t ***segseglist, int segseglistsize, int reqType )
{
		int     i, j, k, l, m ;
        int     setSeg1Forces, setSeg2Forces, coreOnly;
        int     armID12, armID21, armID34, armID43;
        real8   a, MU, NU, Ecore;
        real8   pos1[3], pos2[3], burg[3];
        real8   dx, dy, dz ;
        real8   f1[3], f2[3], f3[3], f4[3];
        real8   sigb[3], extStress[3][3];
		real8	eps;
        Node_t  *node1, *node2, *node3, *node4;
        Param_t *param;

        param    = home->param;
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
 *      Reset all node forces to zero (local and ghost nodes)
 */
        ZeroNodeForces(home, reqType);
		
/*
 *		Code for looping over segseglist or all pair of segments.
 */	
		if (segseglist == NULL) {
			ZeroLongRangeNodeForces(home, FULL);
			
			#ifdef _OPENMP
			#pragma omp parallel for default(none) schedule(dynamic,1) \
				private(i, j , f1 , node1 , armID12 , pos1 , dx , setSeg1Forces) \
				private(k, l , f2 , node2 , armID21 , pos2 , dy , setSeg2Forces) \
				private(       f3 , node3 , armID34 , sigb , dz ) \
				private(burg , f4 , node4 , armID43 , eps  ) \
				shared (a    , MU , Ecore , reqType ) \
				shared (home , NU , param , extStress )
			#endif
			for (i = 0; i < home->newNodeKeyPtr; i++) {

				if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;

				pos1[X] = node1->x;
				pos1[Y] = node1->y;
				pos1[Z] = node1->z;

				for (j = 0; j < node1->numNbrs; j++) {
					node2 = GetNeighborNode(home, node1, j);

					if (node2 == (Node_t *)NULL) {
						printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
						continue;
					}

					// Insures node1 is the node with the lower tag
					if (OrderNodes(node1, node2) >= 0) continue;

					dx = node2->x - pos1[X];
					dy = node2->y - pos1[Y];
					dz = node2->z - pos1[Z];

					ZImage(param, &dx, &dy, &dz);
					
					/* It is possible to have a zero-length segment (created by
					   collision handling).  If we find such a segment, there will
					   be no forces on the segment, so just skip to the next segment. */
		 
					if ((dx*dx + dy*dy + dz*dz) < eps ) continue;

					armID12 = j;
					armID21 = GetArmID(home, node2, node1);
					
					pos2[X] = pos1[X] + dx;
					pos2[Y] = pos1[Y] + dy;
					pos2[Z] = pos1[Z] + dz;

					burg[X] = node1->burgX[armID12];
					burg[Y] = node1->burgY[armID12];
					burg[Z] = node1->burgZ[armID12];

					setSeg1Forces = 1;
					
					/*  If we're doing a partial force calculation, only
						reset forces for this segment if one of the nodes
						is flagged for a force update. */
		 
					if (reqType == PARTIAL) {
						if (((node1->flags & NODE_RESET_FORCES) == 0) &&
							((node2->flags & NODE_RESET_FORCES) == 0)) {
							setSeg1Forces = 0;
						}
					}

					/*  Before calculating the force from other segments,
						calculate forces specific to this segment */

					if (setSeg1Forces) {

						// Add in force due to self stress

						SelfForce(0, MU, NU, burg[X], burg[Y], burg[Z], pos1[X], pos1[Y], pos1[Z],
											 pos2[X], pos2[Y], pos2[Z], a , Ecore, f1  , f2  );
						
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

						// Add in force due to external stress
						
						ExtPKForce(extStress, burg[X], burg[Y], burg[Z],
								   pos1[X], pos1[Y], pos1[Z],
								   pos2[X], pos2[Y], pos2[Z], f1, f2);
						
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

						#ifdef _FEM

							/*  No remote forces are calculated, but if the FEM
								code is hooked in, we need to add in force from
								the FEM stress. */
							sigb[0] = node1->sigbRem[armID12*3];
							sigb[1] = node1->sigbRem[armID12*3+1];
							sigb[2] = node1->sigbRem[armID12*3+2];

							PKForce(sigb, pos1[X], pos1[Y], pos1[Z],
									pos2[X], pos2[Y], pos2[Z], f1, f2);

							AddtoArmForce(node1, armID12, f1);
							AddtoArmForce(node2, armID21, f2);
						#endif
					}

					//  Now compute the force between segment (node1--node2) and all other segments in the simulation.
				
					for (k = 0; k < home->newNodeKeyPtr; k++) {

						if ((node3 = home->nodeKeys[k]) == (Node_t *)NULL) continue;

						for (l = 0; l < node3->numNbrs; l++) {

							node4 = GetNeighborNode(home, node3, l);

							if (node4 == (Node_t *)NULL) {
								printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
								continue;
							}

							// Insures the node with the lower tag is the node3
							if (OrderNodes(node3, node4) >= 0) continue;

							/*  Make sure we don't try to calculate seg/seg forces
							    on a segment with itself, and that we only do
							    forces between a given pair once. */
							if ((node1 == node3) && (node2 == node4)) continue;
							if ( OrderNodes(node3, node1) < 0) continue;
							if ((OrderNodes(node4, node2) < 0) && (node1 == node3)) continue;

							setSeg2Forces = 1;

							/*  If we're doing a partial force calculation, only reset forces for 
								this segment if one of the nodes is flagged for a force update. */

							if (reqType == PARTIAL) {
								if (((node3->flags & NODE_RESET_FORCES) == 0) &&
									((node4->flags & NODE_RESET_FORCES) == 0)) {
									setSeg2Forces = 0;
								}
							}
							
							if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) continue;
							
							armID34 = l;
							armID43 = GetArmID(home, node4, node3);
							
							/*  Calculate the forces between segment (node1--node2)	and segment (node3--node4). */
							ComputeForces(home, node1, node2, node3, node4, f1 , f2 , f3 , f4 );

							if (setSeg1Forces) {
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
							}

							if (setSeg2Forces) {
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
							}
						}
					}
				}
			}
		} else {
			#ifdef _OPENMP
			#pragma omp parallel for default(none) schedule(dynamic,1) \
				private(i ,     f1 , node1 , armID12 , pos1 , dx , setSeg1Forces) \
				private(        f2 , node2 , armID21 , pos2 , dy , setSeg2Forces) \
				private(        f3 , node3 , armID34 , sigb , dz ) \
				private(burg  , f4 , node4 , armID43 ) \
				shared (a     , MU , Ecore , reqType , segseglist , extStress) \
				shared (home  , NU , param , coreOnly, segseglistsize  , eps )
			#endif
			for (i = 0; i < segseglistsize; i++) {        
				node1 = segseglist[i][0];
				node2 = segseglist[i][1];
				node3 = segseglist[i][2];
				node4 = segseglist[i][3];

			
				pos1[X] = node1->x;
				pos1[Y] = node1->y;
				pos1[Z] = node1->z;

				dx = node2->x - pos1[X];
				dy = node2->y - pos1[Y];
				dz = node2->z - pos1[Z];

				ZImage(param, &dx, &dy, &dz);

				/* It is possible to have a zero-length segment (created by
				   collision handling).  If we find such a segment, there will
				   be no forces on the segment, so just skip to the next segment. */
	 
				if ((dx*dx + dy*dy + dz*dz) < eps ) continue;

				armID12 = GetArmID(home,node1,node2);
				armID21 = GetArmID(home,node2,node1);
			
				pos2[X] = pos1[X] + dx;
				pos2[Y] = pos1[Y] + dy;
				pos2[Z] = pos1[Z] + dz;

				setSeg1Forces = 1;

				/*  If we're doing a partial force calculation, only
					reset forces for this segment if one of the nodes
					is flagged for a force update. */
	 
				if (reqType == PARTIAL) {
					if (((node1->flags & NODE_RESET_FORCES) == 0) &&
						((node2->flags & NODE_RESET_FORCES) == 0)) {
						setSeg1Forces = 0;
					}
				}

				if (node1==node3 && node2==node4) {

					burg[X] = node1->burgX[armID12];
					burg[Y] = node1->burgY[armID12];
					burg[Z] = node1->burgZ[armID12];

					/*  Before calculating the force from other segments,
						calculate forces specific to this segment */

					if (setSeg1Forces) {

						// Add in force due to self stress

						SelfForce(coreOnly, MU, NU, burg[X], burg[Y], burg[Z], pos1[X], pos1[Y], pos1[Z],
											 pos2[X], pos2[Y], pos2[Z], a , Ecore, f1  , f2  );
					
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

						// Add in force due to external stress
						ExtPKForce(extStress, burg[X], burg[Y], burg[Z],
								   pos1[X], pos1[Y], pos1[Z],
								   pos2[X], pos2[Y], pos2[Z], f1, f2);
					
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

						#ifdef _FEM

							/*  No remote forces are calculated, but if the FEM
								code is hooked in, we need to add in force from
								the FEM stress. */
							sigb[0] = node1->sigbRem[armID12*3];
							sigb[1] = node1->sigbRem[armID12*3+1];
							sigb[2] = node1->sigbRem[armID12*3+2];

							PKForce(sigb, pos1[X], pos1[Y], pos1[Z],
									pos2[X], pos2[Y], pos2[Z], f1, f2);

							AddtoArmForce(node1, armID12, f1);
							AddtoArmForce(node2, armID21, f2);
						#endif
					}
				} else {

					if (coreOnly) continue;
					setSeg2Forces = 1;

					/*  If we're doing a partial force calculation, only
						reset forces for this segment if one of the nodes
						is flagged for a force update. */

					if (reqType == PARTIAL) {
						if (((node3->flags & NODE_RESET_FORCES) == 0) &&
							((node4->flags & NODE_RESET_FORCES) == 0)) {
							setSeg2Forces = 0;
						}
					}

					if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) continue;
					
					armID34 = GetArmID(home,node3,node4);
					armID43 = GetArmID(home,node4,node3);

					/*  Calculate the forces between segment (node1--node2)	and segment (node3--node4). */
					ComputeForces(home, node1, node2, node3, node4, f1 , f2 , f3 , f4 );
					
					if (setSeg1Forces) {
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
					}

					if (setSeg2Forces) {
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
					}
				}
			}
		}
		
/*
 *      Forces for all segments have been updated, so just
 *      sum up the segment forces to get the nodal forces.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            node1->fX = node1->fxLong;
            node1->fY = node1->fyLong;
            node1->fZ = node1->fzLong;

            for (j = 0; j < node1->numNbrs; j++) {
                node1->fX += node1->armfx[j];
                node1->fY += node1->armfy[j];
                node1->fZ += node1->armfz[j];
            }
        }

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     IncreaseSegSegListSize
 *      Description:   
 *
 *------------------------------------------------------------------------*/
static void IncreaseSegSegListSize(Home_t *home)
{
		int       m, maxSiz;
		Subcyc_t  *subcyc;
		
		subcyc = home->subcyc;
		maxSiz = subcyc->max_SegSegList;
		
		maxSiz +=1000;
		subcyc->ShortRange_SegSegList = realloc (
			subcyc->ShortRange_SegSegList, maxSiz * sizeof(Node_t **));
		
		if ( subcyc->ShortRange_SegSegList==NULL ) {
			puts ("Error (re)allocating memory in IncreaseSegSegListSize function (err 01)");
			exit (1);
		}
		
		for ( m = maxSiz-1000 ; m < maxSiz ; m++) {
			subcyc->ShortRange_SegSegList[m] = malloc(4 * sizeof(Node_t *));
			
			if ( subcyc->ShortRange_SegSegList[m]==NULL ) {
				puts ("Error (re)allocating memory in IncreaseSegSegListSize function (err 02)");
				exit (1);
			}
		}
		
		subcyc->max_SegSegList = maxSiz;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     SegSegList_For_Force
 *      Description:   
 *
 *------------------------------------------------------------------------*/
void SegSegList_For_Force(Home_t *home)
{
		int       i , j , k , l , m ;
		int       cell2Index, nbrCell2Index, nextIndex;
		int       cell2X, cell2Y, cell2Z, cx, cy, cz;
		real8	  dCutOf , dCutOf2;
		real8	  x1, x2, x3, x4;
		real8	  y1, y2, y3, y4;
		real8	  z1, z2, z3, z4;
		real8	  vec1[3], vec2[3];
		Param_t   *param;
		Node_t    *node1 , *node2 , *node3 , *node4 ;
		Subcyc_t  *subcyc;

		param           = home->param;
		subcyc          = home->subcyc;
		dCutOf          = param->cutoff2 ;
		dCutOf2         = dCutOf * dCutOf;
		SortNodes(home , &dCutOf);
		subcyc->Size_SegSegList = 0;
		
		// Allocating memory to ShortRange_SegSegList
		if (subcyc->ShortRange_SegSegList==NULL) {
			subcyc->max_SegSegList = 0.5 * home->newNodeKeyPtr * home->newNodeKeyPtr ;
			
			subcyc->ShortRange_SegSegList = 
				malloc(subcyc->max_SegSegList * sizeof(Node_t **));
				
			if ( subcyc->ShortRange_SegSegList==NULL ) {
				puts ("Error (re)allocating memory in SegSegList_For_Force function (err 01)");
				exit (1);
			}
			
			for ( i = 0 ; i < subcyc->max_SegSegList ; i++) {
				subcyc->ShortRange_SegSegList[i] = malloc(4 * sizeof(Node_t *));
				if ( subcyc->ShortRange_SegSegList[i]==NULL ) {
					puts ("Error (re)allocating memory in SegSegList_For_Force function (err 02)");
					exit (1);
				}
			}
		}
		
        // Start looping through native segments looking for segments pairs.
		for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;

            for (j = 0; j < node1->numNbrs; j++) {

                node2 = GetNeighborNode(home, node1, j);

                if (node2 == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                    continue;
                }
				
				// Ensures the node with the lower tag is the node11
				if (OrderNodes(node1, node2) >= 0) continue;
				
				if (subcyc->Size_SegSegList >= subcyc->max_SegSegList) IncreaseSegSegListSize(home);
				
				subcyc->ShortRange_SegSegList[subcyc->Size_SegSegList][0] = node1;
				subcyc->ShortRange_SegSegList[subcyc->Size_SegSegList][1] = node2;
				subcyc->ShortRange_SegSegList[subcyc->Size_SegSegList][2] = node1;
				subcyc->ShortRange_SegSegList[subcyc->Size_SegSegList][3] = node2;
				subcyc->Size_SegSegList += 1;

				/* Loop through all cell2s neighboring node1.  Only
				   nodes in these neighboring cell2s are candidates. */
				cell2Index = node1->cell2Idx;
				if (cell2Index < 0) continue;
				
				x1 = node1->x;   y1 = node1->y;   z1 = node1->z;
				x2 = node2->x;   y2 = node2->y;   z2 = node2->z;
				
				PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);

				/*   It is possible to have a zero-length segment.
					 If encountered, avoid it.  */
				vec1[0] = x2 - x1;
				vec1[1] = y2 - y1;
				vec1[2] = z2 - z1;
				if (DotProduct(vec1, vec1) < 1.0e-20) continue;

				DecodeCell2Idx(home, cell2Index, &cell2X, &cell2Y, &cell2Z);

				for (cx = cell2X - 1; cx <= cell2X + 1; cx++) {
					for (cy = cell2Y - 1; cy <= cell2Y + 1; cy++) {
						for (cz = cell2Z - 1; cz <= cell2Z + 1; cz++) {
							nbrCell2Index = EncodeCell2Idx(home, cx, cy, cz);

							// Loop though all nodes in the neighbor cell2
							nextIndex = home->cell2[nbrCell2Index];

							while (nextIndex >= 0) {

								node3     = home->cell2QentArray[nextIndex].node;
								nextIndex = home->cell2QentArray[nextIndex].next;

								if (node3 == (Node_t *)NULL) continue;
								
								// Loop over all arms of node3.
								
								for (l = 0; l < node3->numNbrs; l++) {
									node4 = GetNeighborNode(home, node3, l);
									
									// Ensures the node with the lower tag is the node3
									if (OrderNodes(node3, node4) >= 0) continue;
								
									/*  Make sure this isn't a self term (handled above),
										and that we only do this segment pair once.  */
									
									if ((node1 == node3) && (node2 == node4)) continue;
									if ( OrderNodes(node3, node1) < 0) continue;
									if ((OrderNodes(node4, node2) < 0) && (node1 == node3)) continue;
									

									x3 = node3->x;   y3 = node3->y;   z3 = node3->z;
									x4 = node4->x;   y4 = node4->y;   z4 = node4->z;

									PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
									PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

									/*   It is possible to have a zero-length segment.
										 If encountered, avoid it.  */
									vec2[0] = x4 - x3;
									vec2[1] = y4 - y3;
									vec2[2] = z4 - z3;
									if (DotProduct(vec2, vec2) < 1.0e-20) continue;
									
									//  Add information for all nodes if within distance criteria.
									if (subcyc->Size_SegSegList >= subcyc->max_SegSegList) IncreaseSegSegListSize(home);

									subcyc->ShortRange_SegSegList[subcyc->Size_SegSegList][0] = node1;
									subcyc->ShortRange_SegSegList[subcyc->Size_SegSegList][1] = node2;
									subcyc->ShortRange_SegSegList[subcyc->Size_SegSegList][2] = node3;
									subcyc->ShortRange_SegSegList[subcyc->Size_SegSegList][3] = node4;
									subcyc->Size_SegSegList += 1;
									
								} //for (node4)
							} //while (nextnode)
						} //for (cz)
					} //for (cy)
				} //for (cx)
			} // for (node2)
        } //for (node1)
		
}
