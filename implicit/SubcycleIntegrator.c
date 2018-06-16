/*-------------------------------------------------------------------------
 *
 *      Function:     SubcycleIntegrator
 *      Description:  This function time integrates the system one time
 *					  step using a "subcycling" approach. Two groups of
 *					  nodes are established and these two groups are then
 *					  time integrated independently. 
 *
 *-----------------------------------------------------------------------*/

#include "Home.h"
#include "Comm.h"
#include "sys/stat.h"
#include "sys/types.h"

#ifdef _GPU_SUBCYCLE
#include "SubcycleGPU.h"
#endif

#define PI 3.14159265

void SortNodes(Home_t *home, real8 *maxsep);
void ImplicitIntegratorSub0(Home_t *home);
void FlagGroup0Nodes(Home_t *home);
void FlagGroup1Nodes(Home_t *home);
void NodeForceList(Home_t *home, Node_t ***segseglist, int segseglistsize, int reqType );
void TrapezoidIntegratorSub(Home_t *home , int reqType);
void TrapezoidIntegratorSub1(Home_t *home, Node_t ***seglist1, int seglist1size);
void RKFIntegrator(Home_t *home, int reqType);
void SegSegListMaker(Home_t *home, int reqType);
void CommSendCoord(Home_t *home, int subGroup);
void CommSendVelocitySub(Home_t *home, int subGroup);
void ImplicitIntegratorSub(Home_t *home , int reqType);

void AssignSubcycleGroups(Home_t *home, Node_t ***seglistS, int *seglistSsize,
										 int *group1size)
{
		int     eloc, i, j, k, q, arm12, arm21, arm34, arm43;
		int     splitSeg1, splitSeg2, tmpint;
		int     cell2Index, nbrCell2Index, nextIndex;
		int     cell2X, cell2Y, cell2Z, cx, cy, cz;
		int		N, ind, arm, tag1, tag2, tag3, tag4, tagtmp;
		int		flag2, flag3, flag4;
		int		foundmatch, same13, same14, same23, same24;
		int 	hinge;
		real8 	sign1, sign2;
		real8	tmp, dist2, dist2_1, dist2_2, ddist2dt, L1, L2;
		real8	len1, len2;
		real8	x1, x2, x3, x4;
		real8	y1, y2, y3, y4;
		real8	z1, z2, z3, z4;
		real8	xm, ym, zm, xh1, yh1, zh1, xh2, yh2, zh2;
		real8	vXm, vYm, vZm, vXh1, vYh1, vZh1, vXh2, vYh2, vZh2; 
		real8	costheta, theta, sign;
		real8	vec1[3], vec2[3], velvec[3], dr1[3], dr2[3];
		real8	r12[3], r13[3], v12[3], v13[3];
		real8	norm_r12, norm_r13, dot1, dot2, dot3, dot4, dot5;
		real8	D, thetadot;
		real8 	tmp1, tmp2, tmp3, area;

		Node_t  *node1, *node2, *node3, *node4, *nbr, *nodetmp, *node;
		Node_t	*nodem, *nodeh1, *nodeh2;
		Param_t *param;

		param = home->param;

		N = home->newNodeKeyPtr;

		//NEED TO MAKE THIS AN INPUT CONTROL PARAMETER
		//real8	rgroup = 1;
		real8	rgroup = param->rg1;
		real8	thetagroup = 60;
		real8	lgroup = 0.25*param->minSeg;
		real8	area_group = 0.5*param->minSeg*param->minSeg;
		real8 	cell2size = rgroup+param->maxSeg/2;
		//real8	cell2size = 2*param->maxSeg;
		/////////

/*
 *		Initialize all nodes to be in group 0.
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			node->subgroup = 0;
		}
		*group1size = 0;

/*
 *		Add small segments to group 1.
 */
		real8	dx, dy, dz;
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			for (j=0; j<node->numNbrs; j++) {
				nbr = GetNeighborNode(home,node,j);
				if (OrderNodes(node,nbr)>0) continue;
				dx = node->x-nbr->x;
				dy = node->y-nbr->y;
				dz = node->z-nbr->z;
				ZImage(param,&dx,&dy,&dz);
				if (sqrt(dx*dx+dy*dy+dz*dz)<=lgroup) {
					if (node->subgroup==0) {
						node->subgroup = 1;
						(*group1size)++;
					}
					if (nbr->subgroup==0) {
						nbr->subgroup = 1;
						(*group1size)++;
					}
				}
			}
		}

/*
 * 		Sort nodes into cell2s of the grouping radius in size.
 */
		SortNodes(home, &cell2size);

/*
 *      Start looping through native nodes looking for segments pairs...
 */		
		for (tag1 = 0; tag1 < home->newNodeKeyPtr; tag1++) {

            if ((node1 = home->nodeKeys[tag1]) == (Node_t *)NULL) continue;

/*
 *          Loop through all cell2s neighboring the node.  Only
 *          nodes in these neighboring cell2s are candidates for
 *          elastic enhancement.
 */
            cell2Index = node1->cell2Idx;
            if (cell2Index < 0) {
                continue;
            }

            DecodeCell2Idx(home, cell2Index, &cell2X, &cell2Y, &cell2Z);

            for (cx = cell2X - 1; cx <= cell2X + 1; cx++) {
             for (cy = cell2Y - 1; cy <= cell2Y + 1; cy++) {
              for (cz = cell2Z - 1; cz <= cell2Z + 1; cz++) {

                nbrCell2Index = EncodeCell2Idx(home, cx, cy, cz);

/*
 *              Loop though all nodes in the neighbor cell2
 */
                nextIndex = home->cell2[nbrCell2Index];

                while (nextIndex >= 0) {

                    node3 = home->cell2QentArray[nextIndex].node;
                    nextIndex = home->cell2QentArray[nextIndex].next;			

                    if (node3 == (Node_t *)NULL) continue;
/*
 *                  Loop over all arms of node1.
 */
                    for (arm12 = 0; arm12 < node1->numNbrs; arm12++) {

                        node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);

                        if (node2 == (Node_t *)NULL) continue;

                        if ((OrderNodes(node1, node2) >= 0)) {
                            continue;
                        }

/*
 *                      Loop over all arms of node3.
 */
                        for (arm34 = 0; arm34 < node3->numNbrs; arm34++) {

                            node4 = GetNodeFromTag(home, node3->nbrTag[arm34]);

                            if (node4 == (Node_t *)NULL) continue;

/*
 *                      	Ensures the node with the lower tag is the node3
 */
							if (OrderNodes(node3, node4) >= 0) {
								continue;
							}

/*
 *                      	Make sure this isn't a self term (handled above),
 *                      	and that we only do this segment pair once.
 */
							if ((node1 == node3) && (node2 == node4)) {
								continue;
							}

							if ((OrderNodes(node3, node1) < 0)) {
								continue;
							}

							if ((OrderNodes(node4, node2) < 0) &&
								(node1 == node3)) {
								continue;
							}
							if(node1==node4) {
								printf("node1=node4!\n");
							}

                            x1 = node1->x; y1 = node1->y; z1 = node1->z;
                            x2 = node2->x; y2 = node2->y; z2 = node2->z;
                            x3 = node3->x; y3 = node3->y; z3 = node3->z;
                            x4 = node4->x; y4 = node4->y; z4 = node4->z;

                            PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
                            PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
                            PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

/*
 *                          It is possible to have a zero-length segment.
 *                          If encountered, avoid it.
 */
                            vec1[0] = x2 - x1;
                            vec1[1] = y2 - y1;
                            vec1[2] = z2 - z1;
							len1 = sqrt(DotProduct(vec1, vec1));

                            if (len1 < 1.0e-20) {
                                continue;
                            }

                            vec2[0] = x4 - x3;
                            vec2[1] = y4 - y3;
                            vec2[2] = z4 - z3;
							len2 = sqrt(DotProduct(vec2, vec2));

                            if (len2 < 1.0e-20) {
                                continue;
                            }

/*
 *							If two nodes are the same, this is a hinge geometry
 *							and we want to determine the distance between the
 *							free node of the shorter arm and the other segment.
 *							Also store a few things to be used below.
 */
							hinge = 0;
							if ((node1==node3) || 
										(node2==node3) || 
												(node2==node4)) {
								hinge = 1;
								if (node1==node3) {
									sign = 1;
									nodem = node1; nodeh1 = node2; nodeh2 = node4;
									xm = x1; ym = y1; zm = z1;
									xh1 = x2; yh1 = y2; zh1 = z2;
									xh2 = x4; yh2 = y4; zh2 = z4;
									vXm = node1->vX; 	vYm = node1->vY; 	vZm = node1->vZ;
									vXh1 = node2->vX; 	vYh1 = node2->vY; 	vZh1 = node2->vZ;
									vXh2 = node4->vX; 	vYh2 = node4->vY; 	vZh2 = node4->vZ;
								} else if (node2==node3) {
									sign = -1;
									nodem = node2; nodeh1 = node1; nodeh2 = node4;
									xm = x2; ym = y2; zm = z2;
									xh1 = x1; yh1 = y1; zh1 = z1;
									xh2 = x4; yh2 = y4; zh2 = z4;
									vXm = node2->vX; 	vYm = node2->vY; 	vZm = node2->vZ;
									vXh1 = node1->vX; 	vYh1 = node1->vY; 	vZh1 = node1->vZ;
									vXh2 = node4->vX; 	vYh2 = node4->vY; 	vZh2 = node4->vZ;
								} else if (node2==node4) { 
									sign = 1;
									nodem = node2; nodeh1 = node1; nodeh2 = node3;
									xm = x2; ym = y2; zm = z2;
									xh1 = x1; yh1 = y1; zh1 = z1;
									xh2 = x3; yh2 = y3; zh2 = z3;
									vXm = node2->vX; 	vYm = node2->vY; 	vZm = node2->vZ;
									vXh1 = node1->vX; 	vYh1 = node1->vY; 	vZh1 = node1->vZ;
									vXh2 = node3->vX; 	vYh2 = node3->vY; 	vZh2 = node3->vZ;
								}								

								if (len1>len2) {
									GetMinDist2(xm, ym, zm, 0, 0, 0,
												xh1, yh1, zh1, 0, 0, 0,
												xh2, yh2, zh2, 0, 0, 0,
												xh2, yh2, zh2, 0, 0, 0,
												&dist2, &ddist2dt, &L1, &L2); 
								} else {
									GetMinDist2(xm, ym, zm, 0, 0, 0,
												xh2, yh2, zh2, 0, 0, 0,
												xh1, yh1, zh1, 0, 0, 0,
												xh1, yh1, zh1, 0, 0, 0,
												&dist2, &ddist2dt, &L1, &L2); 
								}
							}
/*
 *                          Find the minimum distance between the two segments.
 */
							if (hinge==0) {
								GetMinDist2(x1, y1, z1, 0, 0, 0,
											x2, y2, z2, 0, 0, 0,
											x3, y3, z3, 0, 0, 0,
											x4, y4, z4, 0, 0, 0,
											&dist2, &ddist2dt, &L1, &L2); 
							}

/*
 * 							Test to see if this segment pair belongs to group 1.
 */
							if ( sqrt(dist2) < rgroup ) {
								if (node1->subgroup==0) {
									node1->subgroup = 1;
									(*group1size)++;
								}
								if (node2->subgroup==0) {
									node2->subgroup = 1;
									(*group1size)++;
								}
								if (node3->subgroup==0) {
									node3->subgroup = 1;
									(*group1size)++;
								}
								if (node4->subgroup==0) {
									node4->subgroup = 1;
									(*group1size)++;
								}

								seglistS[*seglistSsize][0] = node1;
								seglistS[*seglistSsize][1] = node2;
								seglistS[*seglistSsize][2] = node3;
								seglistS[*seglistSsize][3] = node4;
								(*seglistSsize)++;
								continue;
							} 

/*
 *							If we have a hinge geometry (segments share a node),
 *							test to see if the angle is smaller than thetagroup
 *							and if the shared node is moving in a direction
 *							that shortens both segments - this is a highly
 *							nonlinear geometry that must be subcycled.
 */

//							if ((node1==node3) || 
//										(node2==node3) || 
//												(node2==node4)) {
							if (hinge) {
								dr1[0] = vec1[0]; 
								dr1[1] = vec1[1]; 
								dr1[2] = vec1[2];
								dr2[0] = vec2[0]; 
								dr2[1] = vec2[1]; 
								dr2[2] = vec2[2];
								NormalizeVec(vec1);
								NormalizeVec(vec2);
								costheta = DotProduct(vec1, vec2);
								if (costheta>1) costheta = 1.0;
								if (costheta<-1) costheta = -1.0;

								r12[0] = xm-xh1; r12[1] = ym-yh1; r12[2] = zm-zh1;
								r13[0] = xm-xh2; r13[1] = ym-yh2; r13[2] = zm-zh2;
								v12[0] = vXm-vXh1; v12[1] = vYm-vYh1; v12[2] = vZm-vZh1;
								v13[0] = vXm-vXh2; v13[1] = vYm-vYh2; v13[2] = vZm-vZh2;
								theta = acos(sign*costheta)*180/PI;

								norm_r12 = sqrt(DotProduct(r12,r12));
								norm_r13 = sqrt(DotProduct(r13,r13));								
								dot1 = DotProduct(r12,r13);
								dot2 = DotProduct(v13,r12);
								dot3 = DotProduct(r13,v12);
								dot4 = DotProduct(v13,r13);
								dot5 = DotProduct(v12,r12);
								D = dot1/norm_r12/norm_r13;
								thetadot = -1/(norm_r12*norm_r13*sqrt(1-D*D))*
									(dot2+dot3-dot1*(dot4/norm_r13/norm_r13
														+dot5/norm_r12/norm_r12));
								if (theta<thetagroup && thetadot>0
														&& nodem->numNbrs==2) {
									if (nodem->subgroup==0) {
										nodem->subgroup = 1;
										(*group1size)++;
									}			
								}
/*								
 *								Also test for the case of "tiny triangles".
 */
								if (Connected(nodeh1,nodeh2,&tmpint)) {
									xvector(dr1[0], dr1[1], dr1[2],
             								dr2[0], dr2[1], dr2[2],
             								&tmp1, &tmp2, &tmp3);
									area = 0.5*sqrt(tmp1*tmp1+tmp2*tmp2+tmp3*tmp3);
									if (area<area_group) {
										if (nodem->subgroup==0) {
											nodem->subgroup = 1;
											(*group1size)++;
										}
										if (nodeh1->subgroup==0) {
											nodeh1->subgroup = 1;
											(*group1size)++;
										}
										if (nodeh2->subgroup==0) {
											nodeh2->subgroup = 1;
											(*group1size)++;
										}
									}
								}										
							}

                        } //for (arm34)
                    } //for (arm12)
                } //while (nextnode)
              } //for (cz)
             } //for (cy)
            } //for (cx)
        } //for (tag1=main node)

/*
 *		Finally, we need to add self terms to the segment lists. All arms that
 *		have at least one node being subcycled need to be added to seglistS. 
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			for (j=0; j<node->numNbrs; j++) {
				nbr = GetNeighborNode(home,node,j);
				if (OrderNodes(node,nbr)>0) continue;
				if (node->subgroup==1 || nbr->subgroup==1) {
					seglistS[*seglistSsize][0] = node;
					seglistS[*seglistSsize][1] = nbr;
					seglistS[*seglistSsize][2] = node;
					seglistS[*seglistSsize][3] = nbr;
					(*seglistSsize)++;
				}
			}
		}

		return;
}

void ForwardProgressCheckNodeB(Home_t *home, int *group1size)
{
		int		i;
		real8	oldV[3], V[3];
		Node_t	*node;
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			if (node->subgroup==0) {
				continue;
			}
			oldV[0] = node->oldvX;
			oldV[1] = node->oldvY;
			oldV[2] = node->oldvZ;
			V[0] = node->vX;
			V[1] = node->vY;
			V[2] = node->vZ;
			if (DotProduct(oldV,V)<0) {
				node->subgroup = 0;
				(*group1size)--;
				//printf("Removed node %i from subcycling due to velocity reversal\n",node->myTag.index);
			}
		}

		return;
}


void SubcycleIntegratorNodeB(Home_t *home)
{
		real8	   totalsubDT, oldDTsub, newDTsub;
		Param_t	  *param;
		Subcyc_t  *subcyc;
		int		   doAll = 1, mobIterError, group1size, cutDT, count;
		int  	   N = home->newNodeKeyPtr, numpairs;
		int		   seglist1size;

		param   = home->param;
		subcyc  = home->subcyc;
/*
 *		Loop over all nodes and assign nodes to group 1 as appropriate. 
 *		Simultaneously assemble two arrays of segment pairs - 1) segment 
 *		pairs that are members of group 1 (both nodes in group 1) and 
 *		also neighbors within rgroup and 2) segment pairs between group 1
 *		and all other segments. 
 */
		int nodecount = 0, i, j;
		Node_t	*node, *nbr;
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			nodecount++;
		}

/*
 *		Allocate memory for the segment pair lists.
 */
		numpairs = (int)ceil(0.5*(real8)N*((real8)N+1)+(real8)N);
		Node_t  ***seglist1 = malloc(numpairs * sizeof(Node_t **));
		real8	**ftotal    = malloc(  N      * sizeof(real8 *));
		real8	**Xold      = malloc(  N      * sizeof(real8 *));
		
		for(i=0;i<numpairs;i++){
			seglist1[i] = malloc(4 * sizeof(Node_t *));
		}
		for(i=0;i<N;i++){
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			
			Xold  [i] = malloc(3 * sizeof(real8));
			ftotal[i] = malloc(3 * sizeof(real8));
			
			Xold[i][0] = node->x;
			Xold[i][1] = node->y;
			Xold[i][2] = node->z;
		}
		seglist1size = 0;

		AssignSubcycleGroups(home, seglist1, &seglist1size, &group1size);

		subcyc->Group1Frac = 100.0 * group1size/ nodecount;
		printf("group 1 size fraction is %e\n", subcyc->Group1Frac);

		// Time integrate group 0 nodes.
		if (group1size<nodecount) {
			ImplicitIntegratorSub0(home);
		}

/*
 *		Calculate the forces for segments not on the subcycle update list.
 */
		FlagGroup1Nodes(home);
        NodeForce(home, PARTIAL);
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			if (node->subgroup==0) continue;
			node->fxLong = 0; node->fyLong = 0; node->fzLong = 0;
			ftotal[i][0] = node->fX;
			ftotal[i][1] = node->fY;
			ftotal[i][2] = node->fZ;
		}

		FlagGroup1Nodes(home);
		NodeForceList(home, seglist1, seglist1size, PARTIAL);
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			if (node->subgroup==0) continue;
			node->fxLong = ftotal[i][0] - node->fX; 
			node->fX     = ftotal[i][0];
			node->fyLong = ftotal[i][1] - node->fY; 
			node->fY     = ftotal[i][1];
			node->fzLong = ftotal[i][2] - node->fZ; 
			node->fZ     = ftotal[i][2];
		}

/*
 *		Time integrate group 1 nodes (subcycle).
 */
		if (home->cycle==0) {
			param->nextDTsub = param->realdt;
		}
		totalsubDT = 0;
		cutDT = 0;
		subcyc->numSubCycle1 = 0;
		if ( group1size > 0 ) {
			if ( group1size == nodecount ) {
				ImplicitIntegrator(home);
				
			} else {
			
				while ( totalsubDT < param->realdt ) {
					if ( totalsubDT + param->nextDTsub > param->realdt ) {
						oldDTsub = param->nextDTsub;
						param->nextDTsub = param->realdt-totalsubDT;
						newDTsub = param->nextDTsub;					
						cutDT = 1;
					}
					TrapezoidIntegratorSub1(home,seglist1,seglist1size);
/*
 *					Test for forward progress.
 */
					ForwardProgressCheckNodeB(home, &group1size);
					if (group1size==0) break;
					//GenerateOutput(home, STAGE_CYCLE);
					totalsubDT = totalsubDT + param->realdtsub;
					//if (param->realdtsub<1e-14) {
					//	param->nextDTsub = 1e-11;					
					//	break;
					//}
					if (cutDT && param->realdtsub==newDTsub) param->nextDTsub = oldDTsub;
					subcyc->numSubCycle1++;
				}
			}
		}
		printf("%i subcycles\n",subcyc->numSubCycle1);
		
/*
 *		Restoring old positions
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
			if (node == (Node_t *)NULL) continue;
			
			node->oldx = Xold[i][0];
			node->oldy = Xold[i][1];
			node->oldz = Xold[i][2];
		}
		
/*
 *		Update the group 0 velocities.
 */
		FlagGroup0Nodes(home);
        NodeForce(home, PARTIAL);
        mobIterError = CalcNodeVelocities(home, 0, !doAll);
        CommSendVelocity(home);

/*
 *		Unflag all nodes.
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			node->flags &= (~NODE_RESET_FORCES);	
		}		

		for(i=0;i<numpairs;i++){
			free(seglist1[i]);
		}
		for(i=0;i<N;i++){
			node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			
			free(ftotal[i]);
			free(Xold  [i]);
		}
		free(seglist1);
		free(ftotal);
		free(Xold  );

}


/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************/

void ForwardProgressCheckForceB(Home_t *home)
{
		int		i;
		real8	oldV[3], V[3];
		Node_t	*node;
		
		TimerStart(home, FORWARD_PROGRES_CHECK);
		
		#ifdef _OPENMP
		#pragma omp parallel for default(none) schedule(static) \
			private(i , node , oldV , V) \
			shared (home  )
		#endif
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
			if (node == (Node_t *)NULL) continue;
			if (node->subgroup == 0) continue;

		//	oldV[0] = node->oldvXsub1;
		//	oldV[1] = node->oldvYsub1;
		//	oldV[2] = node->oldvZsub1;
			oldV[0] = node->oldvX;
			oldV[1] = node->oldvY;
			oldV[2] = node->oldvZ;
			V[0] = node->vX;
			V[1] = node->vY;
			V[2] = node->vZ;
			if (DotProduct(oldV,V)<0) node->subgroup = 0 ;
		}
		
		node = home->ghostNodeQ;
		while (node) {
			if (node->subgroup == 0) {
				node = node->next;
				continue;
			}
			oldV[0] = node->oldvX;
			oldV[1] = node->oldvY;
			oldV[2] = node->oldvZ;
			V[0] = node->vX;
			V[1] = node->vY;
			V[2] = node->vZ;
			if (DotProduct(oldV,V)<0) node->subgroup = 0 ;
			node = node->next;
		}

		TimerStop(home, FORWARD_PROGRES_CHECK);
		return;
}

void ForwardProgressCheckForceB2(Home_t *home, real8 **oldVarr, real8 **oldVghostarr)
{
		int		i;
		real8	oldV[3], V[3];
		Node_t	*node;
		TimerStart(home, FORWARD_PROGRES_CHECK);
		
		#ifdef _OPENMP
		#pragma omp parallel for default(none) schedule(static) \
			private(i , node , oldV , V) \
			shared (home , oldVarr)
		#endif
		for (i=0; i < home->newNodeKeyPtr; i++) {
			node = home->nodeKeys[i];
			if (node == (Node_t *)NULL) continue;
			if (node->subgroup == 0) continue;
			oldV[0] = oldVarr[i][X];
			oldV[1] = oldVarr[i][Y];
			oldV[2] = oldVarr[i][Z];
			V[0] = node->vX;
			V[1] = node->vY;
			V[2] = node->vZ;
			if (DotProduct(oldV,V)<0) node->subgroup = 0 ;
		}
		
		node = home->ghostNodeQ;
		i = 0;
		while (node) {
			if (node->subgroup == 0) {
				node = node->next;
				continue;
			}
			oldV[0] = oldVghostarr[i][X];
			oldV[1] = oldVghostarr[i][Y];
			oldV[2] = oldVghostarr[i][Z];
			V[0] = node->vX;
			V[1] = node->vY;
			V[2] = node->vZ;
			if (DotProduct(oldV,V)<0) node->subgroup = 0 ;
			i++;
			node = node->next;
		}
		
		TimerStop(home, FORWARD_PROGRES_CHECK);
		return;
}

/*******************************************************************************
 *******************************************************************************/
 
void IntegrateSubGroup(Home_t   *home  , int subGroup , 
                       real8  ***Farm0 , int GroupSize)
{
		int       i , j , nSubcyc , cutDT;
		real8	  nextDTsub , totalsubDT, oldDTsub, newDTsub;
		Param_t	  *param;
		Node_t    *node ;
		Subcyc_t  *subcyc;
		
/*
 *		First check if subcycling is needed for this group
 */
		if (GroupSize == 0) return;

		param  = home->param;
		subcyc = home->subcyc;
		
/*
 *		Initializing
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node->subgroup = 1 ;
		}
		node = home->ghostNodeQ;
		while (node) {			
			node->subgroup = 1 ;
			node = node->next;
		}
		
		if      (subGroup == GROUP1) nextDTsub = param->nextDTsub ;
		else if (subGroup == GROUP2) nextDTsub = param->nextDTsub2;
		else if (subGroup == GROUP3) nextDTsub = param->nextDTsub3;
		else if (subGroup == GROUP4) nextDTsub = param->nextDTsub4;
		
		if ( home->cycle==0 || nextDTsub<=0.0 ) nextDTsub = param->realdt;
		totalsubDT = 0;
		cutDT      = 0;
		nSubcyc    = 0;
				
/*
 *		Subcycle until the subcycle group time (totalsubDT) catches up
 *		with the global group time (realdt). Note that nodal forces 
 *		will reset to zero when subcycling is performed
 */
		while ( totalsubDT < param->realdt ) {
			if ( totalsubDT + nextDTsub > param->realdt ) {
				oldDTsub  = nextDTsub;
				nextDTsub = param->realdt - totalsubDT;
				newDTsub  = nextDTsub;					
				cutDT     = 1;
				
				if      (subGroup == GROUP1) param->nextDTsub  = nextDTsub;
				else if (subGroup == GROUP2) param->nextDTsub2 = nextDTsub;
				else if (subGroup == GROUP3) param->nextDTsub3 = nextDTsub;
				else if (subGroup == GROUP4) param->nextDTsub4 = nextDTsub;
			}
			
			if        ((strcmp(param->subInteg0Integ1, "Exp-Exp") == 0) || 
					   (strcmp(param->subInteg0Integ1, "RKF-Exp") == 0) ){
				TrapezoidIntegratorSub(home, subGroup);
			} else if ((strcmp(param->subInteg0Integ1, "Exp-RKF") == 0) || 
					   (strcmp(param->subInteg0Integ1, "RKF-RKF") == 0) ){
				RKFIntegrator(home, subGroup);
			}
			
			//if (nSubcyc > 3) ForwardProgressCheckForceB(home);
			nSubcyc++;

			if        (subGroup == GROUP1) {
				if (cutDT && param->realdtsub == newDTsub) param->nextDTsub = oldDTsub;
				nextDTsub    = param->nextDTsub ;
				totalsubDT   = totalsubDT + param->realdtsub;
				subcyc->numSubCycle1 = nSubcyc;

			} else if (subGroup == GROUP2) {
				if (cutDT && param->realdtsub2 == newDTsub) param->nextDTsub2 = oldDTsub;
				nextDTsub    = param->nextDTsub2 ;
				totalsubDT   = totalsubDT + param->realdtsub2;
				subcyc->numSubCycle2 = nSubcyc;
				
			} else if (subGroup == GROUP3) {
				if (cutDT && param->realdtsub3 == newDTsub) param->nextDTsub3 = oldDTsub;
				nextDTsub    = param->nextDTsub3 ;
				totalsubDT   = totalsubDT + param->realdtsub3;
				subcyc->numSubCycle3 = nSubcyc;
				
			} else if (subGroup == GROUP4) {
				if (cutDT && param->realdtsub4 == newDTsub) param->nextDTsub4 = oldDTsub;
				nextDTsub    = param->nextDTsub4 ;
				totalsubDT   = totalsubDT + param->realdtsub4;
				subcyc->numSubCycle4 = nSubcyc;
				
			}
				
		//	if (nSubcyc % 5 ==0) printf("Group= %d , nSubcyc= %d\n" , subGroup , nSubcyc);
		}
			
/*
 *		Summing up nodal and arm forces of the groups
 */		
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			node->olderfX += node->fX;
			node->olderfY += node->fY;
			node->olderfZ += node->fZ;
			
			for (j = 0 ; j < node->numNbrs ; j++) {
				Farm0[i][j][0] += node->armfx[j];
				Farm0[i][j][1] += node->armfy[j];
				Farm0[i][j][2] += node->armfz[j];
			}
		}
		node = home->ghostNodeQ;
		while (node) {			
			node->olderfX += node->fX;
			node->olderfY += node->fY;
			node->olderfZ += node->fZ;
			node = node->next;
		}
}

/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************/

void FlagSubcycleNodes(Home_t *home, int subGroup)
{
		int         i, j ,nSeg, nSegSeg;
        Segm_t      *SegList;
        SegSeg_t    *SegSegList;
		Node_t      *node;
		Subcyc_t    *subcyc;
		
		TimerStart(home, FLAG_SUBCYC_NODES);
		subcyc = home->subcyc;

/*
 *      In the FULL case, just flag everything and return.
 */
        if (subGroup == FULL) {
		    #ifdef _OPENMP
			#pragma omp parallel for default(none) schedule(static) \
				private(i , j , node) \
				shared (home  )
			#endif
			for (i=0; i < home->newNodeKeyPtr; i++) {
			    if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
                node->subgroup = 1;
				
                for (j = 0 ; j < 5 ; j++) node->CommSend[j] = 1;
		    }

		    node = home->ghostNodeQ;
		    while (node) {
			    node->subgroup = 1;
			    node = node->next;
		    }
            return;
        }
/*
 *      Initially set all native nodes for subcycling.
 *		Ghost nodes are initially unflaged and then flaged below if appropriate
 */
		
		#ifdef _OPENMP
		#pragma omp parallel for default(none) schedule(static) \
			private(i , node) \
			shared (home  )
		#endif
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            node->subgroup = 1;
        //    node->CommSend = 0;
		}

		node = home->ghostNodeQ;
		while (node) {
			node->subgroup = 0;
			node = node->next;
		}

/*
 *      Loop over the appropriate seg and seg-seg lists and flag all
 *      nodes for subcycling.
 */
		if (subGroup == GROUP0) {
			SegList    = subcyc->SegListG0;
			nSeg       = subcyc->SegListG0_cnt;
			SegSegList = subcyc->SegSegListG0;
			nSegSeg    = subcyc->SegSegListG0_cnt;
			
		} else if (subGroup == GROUP1) {
			SegList    = subcyc->SegListG1;
			nSeg       = subcyc->SegListG1_cnt;
			SegSegList = subcyc->SegSegListG1;
			nSegSeg    = subcyc->SegSegListG1_cnt;
			
		} else if (subGroup == GROUP2) {
			nSeg       = 0;
			SegSegList = subcyc->SegSegListG2;
			nSegSeg    = subcyc->SegSegListG2_cnt;
			
		} else if (subGroup == GROUP3) {
			nSeg       = 0;
			SegSegList = subcyc->SegSegListG3;
			nSegSeg    = subcyc->SegSegListG3_cnt;
			
		} else if (subGroup == GROUP4) {
			nSeg       = 0;
			SegSegList = subcyc->SegSegListG4;
			nSegSeg    = subcyc->SegSegListG4_cnt;
		}

		#ifdef _OPENMP
		#pragma omp parallel for default(none) schedule(dynamic,1) \
			private(i) \
			shared (SegList , nSeg)
		#endif
		for (i=0; i<nSeg; i++) {
			if (SegList[i].flag == 1) {
				SegList[i].seg->node1->subgroup = 1;
				SegList[i].seg->node2->subgroup = 1;
			}
		}

		#ifdef _OPENMP
		#pragma omp parallel for default(none) schedule(dynamic,1) \
			private(i) \
			shared (SegSegList , nSegSeg)
		#endif
		for (i=0; i<nSegSeg; i++) {
			if (SegSegList[i].flag == 1) { 
				SegSegList[i].seg1->node1->subgroup = 1;
				SegSegList[i].seg1->node2->subgroup = 1;
				SegSegList[i].seg2->node1->subgroup = 1;
				SegSegList[i].seg2->node2->subgroup = 1;
			}
		}
	
/*
 *      Loop over the ghost seg list and flag nodes on segments
 *      in the appropriate subGroup for communication.
 */
	/*	for (i=0; i<SegListGhost_cnt; i++) {
			if (SegListGhost[i].flag == subGroup) { 
				SegListGhost[i].seg->node1->CommSend = 1;
				SegListGhost[i].seg->node2->CommSend = 1;
				SegListGhost[i].seg->node1->subgroup = 1;
				SegListGhost[i].seg->node2->subgroup = 1;
			}
		}*/

		TimerStop(home, FLAG_SUBCYC_NODES);
}

/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************/

void SubcycleIntegratorForceB(Home_t *home)
{
		int       TotalSize , doAll = 1 , mobIterError;
		int       GroupSizeLocal[5] , GroupSizeGlobal[5];
		int       i , j;
		Subcyc_t  *subcyc;
		Param_t	  *param;
		Node_t    *node ;
		struct    timeval WCT00 , WCT0 , WCT1 , WCT2 , WCT3 , WCT4;
		struct    timeval WCTold, WCTnew;
		real8     dWCT0 , dWCT1, dWCT2, dWCT3, dWCT4;

		param  = home->param;
		subcyc = home->subcyc;
		
/*
 *      Store the current positions.
 */
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node->oldervX = node->vX;   node->olderx = node->x;
			node->oldervY = node->vY;   node->oldery = node->y;
			node->oldervZ = node->vZ;   node->olderz = node->z;
		}
		node = home->ghostNodeQ;
		while (node) {
			node->oldervX = node->vX;   node->olderx = node->x;
			node->oldervY = node->vY;   node->oldery = node->y;
			node->oldervZ = node->vZ;   node->olderz = node->z;
			
			node = node->next;
		}
		

/*
 *		Making group0 through group4 lists. These lists contain segment
 *		and segment pair information used for force calculations. 
 */
		SegSegListMaker(home, FULL);
		GroupSizeLocal[0] = subcyc->SegSegListG0_cnt + subcyc->SegListG0_cnt;
		GroupSizeLocal[1] = subcyc->SegSegListG1_cnt + subcyc->SegListG1_cnt;
		GroupSizeLocal[2] = subcyc->SegSegListG2_cnt ;
		GroupSizeLocal[3] = subcyc->SegSegListG3_cnt ;
		GroupSizeLocal[4] = subcyc->SegSegListG4_cnt ;
#if PARALLEL		
		MPI_Allreduce(GroupSizeLocal, GroupSizeGlobal, 5, 
		              MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
		for (i=0 ; i<5 ; i++) GroupSizeGlobal[i] = GroupSizeLocal[i];
#endif		
		TotalSize  = GroupSizeGlobal[0] + GroupSizeGlobal[1] + GroupSizeGlobal[2]
		           + GroupSizeGlobal[3] + GroupSizeGlobal[4] ;
		subcyc->Group1Frac = 1.0 * GroupSizeGlobal[1] / TotalSize;
		subcyc->Group2Frac = 1.0 * GroupSizeGlobal[2] / TotalSize;
		subcyc->Group3Frac = 1.0 * GroupSizeGlobal[3] / TotalSize;
		subcyc->Group4Frac = 1.0 * GroupSizeGlobal[4] / TotalSize;
		
	//	if (home->myDomain == 0) gettimeofday( &WCT00, NULL );



#ifdef _GPU_SUBCYCLE
/*
 * 		GPU SUBCYCLING: Call the GPU subcycling integrator if this
 *      option has been selected in the control file. Note that there
 *      are currently some limitations to the use of the GPU subcycling,
 *      e.g. it can only be used in serial mode. See SubcycleGPU.cpp for
 *      more information.
 * 
 * 		Call the GPU subcycling integrator if selected, otherwise just
 *      skip this part and perform the integration using the regular 
 *      CPU version.
 */	
		if (strcmp(param->subInteg0Integ1, "GPU") == 0) {
/*
 *			Perform subcycling on the GPU
 */		
			TimerStart(home, SUBCYCLING_GPU);
			SubcycleIntegratorGPU(home);
			TimerStop(home, SUBCYCLING_GPU);
		
/*
 *			We are done with subcycling. Now restore the old positions
 *			needed for topological changes.
 */		
			CommSendCoord(home, FULL);
			
			for (i=0; i < home->newNodeKeyPtr; i++) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				node->oldx = node->olderx;   node->oldvX = node->oldervX;
				node->oldy = node->oldery;   node->oldvY = node->oldervY;
				node->oldz = node->olderz;   node->oldvZ = node->oldervZ;
			}
			
			return;
		}
		
#else
		if (strcmp(param->subInteg0Integ1, "GPU") == 0) {
			Fatal("GPU subcycling can only be used with compile flag -D_GPU_SUBCYCLE");
		}
#endif


/*
 *		Flag the nodes necessary for GROUP0 subcycling and update the forces
 *		and velocities.
 */
		FlagSubcycleNodes(home, GROUP0);

        NodeForce(home, GROUP0);
        mobIterError = CalcNodeVelocities(home, 0, doAll);
        CommSendVelocitySub(home, GROUP0);
		
	//	printf("myDomain = %d , GroupSizeLocal = %7d %5d %5d %5d %5d\n", home->myDomain, GroupSizeLocal[0],
	//	        GroupSizeLocal[1], GroupSizeLocal[2], GroupSizeLocal[3], GroupSizeLocal[4]);

/*
 *		Time integrate group0 forces.
 */
		if        ((strcmp(param->subInteg0Integ1, "Imp-Exp") == 0) || 
				   (strcmp(param->subInteg0Integ1, "Imp-Imp") == 0) ){
			if (GroupSizeGlobal[0] > 0) ImplicitIntegratorSub(home , GROUP0);
			else                        ImplicitIntegrator   (home );
			
		} else if ((strcmp(param->subInteg0Integ1, "Exp-Exp") == 0) || 
				   (strcmp(param->subInteg0Integ1, "Exp-Imp") == 0) ){
			if (GroupSizeGlobal[0] > 0) TrapezoidIntegratorSub(home, GROUP0);
			else                        TrapezoidIntegrator   (home );
		
		} else if ((strcmp(param->subInteg0Integ1, "RKF-Exp") == 0) || 
				   (strcmp(param->subInteg0Integ1, "RKF-Imp") == 0) || 
				   (strcmp(param->subInteg0Integ1, "RKF-RKF") == 0) ){
			RKFIntegrator(home, GROUP0);
			
		}
		
/*
 *		It is likely that the size of group 4 changed during the 
 *		time-integration of group 0. So adjust its global value 
 *		before proceeding with subcycling.
 */
		if ((param->nTry > 0) && (GroupSizeGlobal[0] > 0)) {
#if PARALLEL
			MPI_Allreduce(&(subcyc->SegSegListG4_cnt), 
			              &(GroupSizeGlobal[4]), 1, 
		                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
			GroupSizeGlobal[4] = subcyc->SegSegListG4_cnt;
#endif
			subcyc->Group4Frac = 1.0 * GroupSizeGlobal[4] / TotalSize;
		}
		
/*
 *		Storing nodal forces corresponding to group0. 
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node->olderfX = node->fX;
			node->olderfY = node->fY;
			node->olderfZ = node->fZ;
		}
		node = home->ghostNodeQ;
		while (node) {			
			node->olderfX = node->fX;
			node->olderfY = node->fY;
			node->olderfZ = node->fZ;
			node = node->next;
		}
		

#ifdef _PRECOMPUTE_SUB_FMM
/*
 * 		Pre-compute the FMM contribution on each segment of group 1
 * 		before entering the subcycling loop. This can significantly
 * 		reduce the computation time of group 1.
 */
		if (param->forceCutOff == 0 && param->fmEnabled) {
			subcyc->sigbFMM = (real8**)malloc(subcyc->SegListG1_cnt*sizeof(real8*));
			RemoteSigbSub(home);
		}
#endif
		

/*
 *		Time integrate group 1, 2, 3 and 4 interactions (subcycle).
 */		
		
	//	if (home->myDomain == 0) {
	//		gettimeofday( &WCT0, NULL );
	//		dWCT0 = (WCT0.tv_sec - WCT00.tv_sec) + (WCT0.tv_usec - WCT00.tv_usec)/1000000.0;
	//	}
/////////////////////////////////////////////////////////////////////////////////////////////
//      NEW SUBCYCLE CODE
        real8   subTime1, subTime2, subTime3, subTime4;
        real8   totalsubDT, nextDTsub, oldDTsub, newDTsub;
        int     subGroup, cutDT;

        //Initialize the wall clock times to zero for subgroups
    /*    dWCT1 = 0.0; 
        dWCT2 = 0.0; 
        dWCT3 = 0.0; 
        dWCT4 = 0.0;*/

        //Initialize the time for each group based on whether it has any forces in it
        if (GroupSizeGlobal[1] > 0) subTime1 = 0.0;
        else                        subTime1 = param->realdt;
		
        if (GroupSizeGlobal[2] > 0) subTime2 = 0.0;
        else                        subTime2 = param->realdt;
		
        if (GroupSizeGlobal[3] > 0) subTime3 = 0.0;
        else                        subTime3 = param->realdt;
		
        if (GroupSizeGlobal[4] > 0) subTime4 = 0.0;
        else                        subTime4 = param->realdt;
		
        //Initialize some other stuff
		if ( home->cycle==0 ) nextDTsub = param->realdt;
		totalsubDT = 0.0;
        subcyc->numSubCycle1 = 0;
        subcyc->numSubCycle2 = 0;
        subcyc->numSubCycle3 = 0;
        subcyc->numSubCycle4 = 0;
        int oldGroup = -1;
        int nSubcyc;
        
        //Allocate memory for array where "old" velocities can be stored,
        //to be used with the forward progress check.
        real8 **oldV, **oldVghost;
        //Native nodes
		oldV = malloc(home->newNodeKeyPtr * sizeof(real8 *));
		for ( i = 0 ; i < home->newNodeKeyPtr ; i++ ) {
			oldV[i] = malloc(3 * sizeof(real8));
		}
        //Ghost nodes
        int ghostsize1 = home->newNodeKeyPtr, ghostsize2;
        if (ghostsize1 == 0) ghostsize1 = 1;
		oldVghost = malloc(ghostsize1 * sizeof(real8 *));
        i = 0;
        node = home->ghostNodeQ;
        while (node) {
            oldVghost[i] = malloc(3 * sizeof(real8));
            i++;
            if (i >= ghostsize1) {
                ghostsize1 *= 2;			
                oldVghost = realloc(oldVghost, ghostsize1 * sizeof(real8 *));
            }
	        node = node->next;
        }
        ghostsize2 = i-1;		
		
/*
 *		Subcycle until the subcycle group times (subTimei) catch up
 *		with the global group time (realdt). Note that nodal forces 
 *		will reset to zero when subcycling is performed
 */
		while ( subTime1 < param->realdt || subTime2 < param->realdt ||
                subTime3 < param->realdt || subTime4 < param->realdt ) {
            cutDT = 0;

            //The group that is furthest behind goes first
            if        ( subTime4 <= subTime3 && subTime4 <= subTime2 && subTime4 <= subTime1 ) {
                subGroup   = GROUP4;
                nextDTsub  = param->nextDTsub4;
                totalsubDT = subTime4;
				
            } else if ( subTime3 < subTime4 && subTime3 <= subTime2 && subTime3 <= subTime1 ) {
                subGroup   = GROUP3;
                nextDTsub  = param->nextDTsub3;
                totalsubDT = subTime3;
				
            } else if ( subTime2 < subTime4 && subTime2 < subTime3 && subTime2 <= subTime1 ) {
                subGroup   = GROUP2;
                nextDTsub  = param->nextDTsub2;
                totalsubDT = subTime2;
				
            } else {
                subGroup   = GROUP1;
                nextDTsub  = param->nextDTsub;
                totalsubDT = subTime1;
            }
    
            //If we switched groups, reset subcycle count and toggle all nodes on
            //for subcycling. Also do bookkeeping with the wall clock times.
            if (subGroup != oldGroup) {

                nSubcyc = 0;
                
				//Flag appropriate nodes for subcycling, they may be unflagged during the forward
				//progress check.
				FlagSubcycleNodes(home, subGroup);

                //Update positions of ghost nodes so they are correct for the new subcycle group.
                CommSendCoord(home, subGroup);

                //Update the forces and velocities using the new group
		        NodeForce(home, subGroup);
		        mobIterError = CalcNodeVelocities(home, 0, doAll);
                CommSendVelocitySub(home, subGroup);

			/*    if (home->myDomain == 0) {
                    gettimeofday( &WCTnew, NULL );

                    if        (oldGroup == GROUP1) {
				        dWCT1 += (WCTnew.tv_sec - WCTold.tv_sec) + (WCTnew.tv_usec - WCTold.tv_usec)/1000000.0;
			        } else if (oldGroup == GROUP2) {
				        dWCT2 += (WCTnew.tv_sec - WCTold.tv_sec) + (WCTnew.tv_usec - WCTold.tv_usec)/1000000.0;
			        } else if (oldGroup == GROUP3) {
				        dWCT3 += (WCTnew.tv_sec - WCTold.tv_sec) + (WCTnew.tv_usec - WCTold.tv_usec)/1000000.0;
			        } else if (oldGroup == GROUP4) {
				        dWCT4 += (WCTnew.tv_sec - WCTold.tv_sec) + (WCTnew.tv_usec - WCTold.tv_usec)/1000000.0;
			        }

                    WCTold = WCTnew;
                }*/
            }

            oldGroup = subGroup;

            //Make sure we don't pass the global group in time
			if ( totalsubDT + nextDTsub > param->realdt ) {
				oldDTsub  = nextDTsub;
				nextDTsub = param->realdt - totalsubDT;
				newDTsub  = nextDTsub;					
				cutDT     = 1;
				
				if      (subGroup == GROUP1) param->nextDTsub  = nextDTsub;
				else if (subGroup == GROUP2) param->nextDTsub2 = nextDTsub;
				else if (subGroup == GROUP3) param->nextDTsub3 = nextDTsub;
				else if (subGroup == GROUP4) param->nextDTsub4 = nextDTsub;
			}
			
            //Time integrate the chosen group for one subcycle
			if        ((strcmp(param->subInteg0Integ1, "Exp-Exp") == 0) || 
					   (strcmp(param->subInteg0Integ1, "RKF-Exp") == 0) ){
				TrapezoidIntegratorSub(home, subGroup);
			} else if ((strcmp(param->subInteg0Integ1, "Exp-RKF") == 0) || 
					   (strcmp(param->subInteg0Integ1, "RKF-RKF") == 0) ){
				RKFIntegrator(home, subGroup);
			}
			
            //Perform forward progress check
		//	if (nSubcyc > 3) ForwardProgressCheckForceB2(home,oldV,oldVghost);
			if (nSubcyc > 3) ForwardProgressCheckForceB(home);
			nSubcyc++;

            //Do bookkeeping on the time step and number of subcycles
			if        (subGroup == GROUP1) {
				if (cutDT && param->realdtsub == newDTsub) param->nextDTsub = oldDTsub;
				subTime1 += param->realdtsub;
				subcyc->numSubCycle1++;
			} else if (subGroup == GROUP2) {
				if (cutDT && param->realdtsub2 == newDTsub) param->nextDTsub2 = oldDTsub;
				subTime2 += param->realdtsub2;
				subcyc->numSubCycle2++;
			} else if (subGroup == GROUP3) {
				if (cutDT && param->realdtsub3 == newDTsub) param->nextDTsub3 = oldDTsub;
				subTime3 += param->realdtsub3;
				subcyc->numSubCycle3++;
			} else if (subGroup == GROUP4) {
				if (cutDT && param->realdtsub4 == newDTsub) param->nextDTsub4 = oldDTsub;
				subTime4 += param->realdtsub4;
				subcyc->numSubCycle4++;
			}
		}

///////////////////////////////////////////////////////////////////////////////////////////

/*
//      Old subcycle code which caused "The Jumble Bug"	
		if ( (GroupSizeGlobal[0] > 0) && (GroupSizeGlobal[0] < TotalSize) ) {
			IntegrateSubGroup(home , GROUP1 , Farm0 , GroupSizeGlobal[1]); 
			if (home->myDomain == 0) gettimeofday( &WCT1, NULL );
		
			IntegrateSubGroup(home , GROUP2 , Farm0 , GroupSizeGlobal[2]);
			if (home->myDomain == 0) gettimeofday( &WCT2, NULL );
			
			IntegrateSubGroup(home , GROUP3 , Farm0 , GroupSizeGlobal[3]);
			if (home->myDomain == 0) gettimeofday( &WCT3, NULL );
			
			IntegrateSubGroup(home , GROUP4 , Farm0 , GroupSizeGlobal[4]);
			if (home->myDomain == 0) {
				gettimeofday( &WCT4, NULL );
				dWCT1 = (WCT1.tv_sec - WCT0.tv_sec) + (WCT1.tv_usec - WCT0.tv_usec)/1000000.0;
				dWCT2 = (WCT2.tv_sec - WCT1.tv_sec) + (WCT2.tv_usec - WCT1.tv_usec)/1000000.0;
				dWCT3 = (WCT3.tv_sec - WCT2.tv_sec) + (WCT3.tv_usec - WCT2.tv_usec)/1000000.0;
				dWCT4 = (WCT4.tv_sec - WCT3.tv_sec) + (WCT4.tv_usec - WCT3.tv_usec)/1000000.0;
			}			
*/

/*
 *      Make sure all of the ghost node positions are up to date.
 */
       // FlagSubcycleNodes(home, FULL);
        CommSendCoord(home, FULL);
			
/*
 *			We are done with subcycling. Now restore the old positions
 *			and update the mobilities
 */		
			for (i=0; i < home->newNodeKeyPtr; i++) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				node->oldx = node->olderx;   node->oldvX = node->oldervX;
				node->oldy = node->oldery;   node->oldvY = node->oldervY;
				node->oldz = node->olderz;   node->oldvZ = node->oldervZ;
				
			/*	node->fX   = node->olderfX;
				node->fY   = node->olderfY;
				node->fZ   = node->olderfZ;
			*/
			}
			
			node = home->ghostNodeQ;
			while (node) {
				node->oldx = node->olderx;   node->oldvX = node->oldervX;
				node->oldy = node->oldery;   node->oldvY = node->oldervY;
				node->oldz = node->olderz;   node->oldvZ = node->oldervZ;
				
			/*	node->fX   = node->olderfX;
				node->fY   = node->olderfY;
				node->fZ   = node->olderfZ;*/
				node = node->next;
			}
			
			NodeForce(home, FULL);
			mobIterError = CalcNodeVelocities(home, 0, doAll);
			CommSendVelocity(home);
//		}
		
	/*	if (home->myDomain == 0) {
		//	printf("%i %i %i %i subcycles\n",subcyc->numSubCycle1, subcyc->numSubCycle2, subcyc->numSubCycle3, subcyc->numSubCycle4);
			FILE   *fp;
			fp = fopen("dWCT", "a");
			fprintf(fp, "%f , %e %i %f , %e %i %f , %e %i %f  , %e %i %f\n", dWCT0 , 
						subcyc->Group1Frac , subcyc->numSubCycle1 , dWCT1 ,
						subcyc->Group2Frac , subcyc->numSubCycle2 , dWCT2 ,
						subcyc->Group3Frac , subcyc->numSubCycle3 , dWCT3 ,
						subcyc->Group4Frac , subcyc->numSubCycle4 , dWCT4 );
			
			fp = fopen("subcycling.txt", "a");
			fprintf(fp, "%d %d %d %d %d %d\n", home->cycle, 1,
					subcyc->numSubCycle1, subcyc->numSubCycle2, 
					subcyc->numSubCycle3, subcyc->numSubCycle4);
			fclose(fp);
		} */
		
/*
 *		Freeing memory
 */ 
#ifdef _PRECOMPUTE_SUB_FMM		
		if (param->forceCutOff == 0 && param->fmEnabled) {
			for (i = 0; i < subcyc->SegListG1_cnt; i++) {
				free(subcyc->sigbFMM[i]);
			}
			free(subcyc->sigbFMM);
		}
#endif
		
		for (i=0; i < home->newNodeKeyPtr; i++) {
			free( oldV[i] );
		}
		free( oldV );

        for (i=0; i < ghostsize2; i++) {
			free( oldVghost[i] );
		}
		free( oldVghost );

}

