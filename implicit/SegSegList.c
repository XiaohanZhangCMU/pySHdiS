#include "Home.h"

int CellPriority(Home_t *home, int cellID1, int cellID2);

void GetMinDistSegSeg(Home_t *home  , Node_t *node1 , Node_t *node2 , 
                      Node_t *node3 , Node_t *node4 , real8  *dist2 )
{
		real8     x1 , x2 , x3 , x4 , xm , xh1 , xh2 , dx21 , dx43 , len21s , L1;
		real8     y1 , y2 , y3 , y4 , ym , yh1 , yh2 , dy21 , dy43 , len43s , L2;
		real8     z1 , z2 , z3 , z4 , zm , zh1 , zh2 , dz21 , dz43 , ddist2dt   ;
		Param_t  *param;
		
		param = home->param;
		
#ifdef _GPU_SUBCYCLE
/*				
 * 		Initially, assign all segments / segments interactions 
 * 		to group 0 when using the GPU subcycle integrator to 
 * 		maximize code performance. The actual group assignement
 * 		will be directly performed on the GPU (much faster).
 */
		*dist2 = 2.0 * param->rg4 * param->rg4;
		return;
#endif
		
		*dist2 = -1;
		if (!param->elasticinteraction) return;
		
		x1 = node1->x; y1 = node1->y; z1 = node1->z;
		x2 = node2->x; y2 = node2->y; z2 = node2->z;
		x3 = node3->x; y3 = node3->y; z3 = node3->z;
		x4 = node4->x; y4 = node4->y; z4 = node4->z;

		PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
		PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
		PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

/*
 *		It is possible to have a zero-length segment.
 *		If encountered, avoid it.
 */
		dx21   = x2 - x1;
		dy21   = y2 - y1;
		dz21   = z2 - z1;
		len21s = dx21*dx21 + dy21*dy21 + dz21*dz21;
		if (len21s < 1.0e-20) return;

		dx43 = x4 - x3;
		dy43 = y4 - y3;
		dz43 = z4 - z3;
		len43s = dx43*dx43 + dy43*dy43 + dz43*dz43;
		if (len43s < 1.0e-20) return;

/*
 *		Find the minimum distance between the two segments.
 *
 *		If two nodes are the same, this is a hinge geometry
 *		and we want to determine the distance between the
 *		free node of the shorter arm and the other segment.
 *		Also store a few things to be used below.
 */	
		int hinge = 0;
		if        ((node1->myTag.domainID == node3->myTag.domainID) &&
		           (node1->myTag.index    == node3->myTag.index   )) {
			hinge = 1;
			
		} else if ((node2->myTag.domainID == node3->myTag.domainID) &&
		           (node2->myTag.index    == node3->myTag.index   )) {
			hinge = 2;
			
		} else if ((node2->myTag.domainID == node4->myTag.domainID) &&
		           (node2->myTag.index    == node4->myTag.index   )) {
			hinge = 3;
			
		} else if ((node1->myTag.domainID == node4->myTag.domainID) &&
		           (node1->myTag.index    == node4->myTag.index   )) {
			hinge = 4;
			
		}
		
		if (!hinge) {
			GetMinDist2(x1, y1, z1, 0, 0, 0, x2, y2, z2, 0, 0, 0, 
			            x3, y3, z3, 0, 0, 0, x4, y4, z4, 0, 0, 0, 
			            dist2, &ddist2dt, &L1, &L2);
		} else {
			if (hinge == 1) {
				xm  = x1 ;  ym  = y1 ;  zm  = z1;
				xh1 = x2 ;  yh1 = y2 ;  zh1 = z2;
				xh2 = x4 ;  yh2 = y4 ;  zh2 = z4;
				
			} else if (hinge == 2) {
				xm  = x2 ;  ym  = y2 ;  zm  = z2;
				xh1 = x1 ;  yh1 = y1 ;  zh1 = z1;
				xh2 = x4 ;  yh2 = y4 ;  zh2 = z4;
				
			} else if (hinge == 3) {
				xm  = x2 ;  ym  = y2 ;  zm  = z2;
				xh1 = x1 ;  yh1 = y1 ;  zh1 = z1;
				xh2 = x3 ;  yh2 = y3 ;  zh2 = z3;
				
			} else if (hinge == 4) {
				xm  = x1 ;  ym  = y1 ;  zm  = z1;
				xh1 = x2 ;  yh1 = y2 ;  zh1 = z2;
				xh2 = x3 ;  yh2 = y3 ;  zh2 = z3;
			}
			
			if (len21s>len43s) {
				GetMinDist2(xm , ym , zm , 0, 0, 0, xh1, yh1, zh1, 0, 0, 0,
				            xh2, yh2, zh2, 0, 0, 0, xh2, yh2, zh2, 0, 0, 0,
				            dist2, &ddist2dt, &L1, &L2); 
							
			} else {
				GetMinDist2(xm , ym , zm , 0, 0, 0, xh2, yh2, zh2, 0, 0, 0, 
				            xh1, yh1, zh1, 0, 0, 0, xh1, yh1, zh1, 0, 0, 0, 
				            dist2, &ddist2dt, &L1, &L2); 
			}
		}
		
		if ((*dist2 > param->cutoff2 * param->cutoff2) &&
		    (param->forceCutOff)) *dist2 = -1.0;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void InitSegLists(Home_t *home)
{
	Subcyc_t  *subcyc;
	subcyc  =  home->subcyc;

	if ( subcyc->SegListG0 != (Segm_t *)NULL ) free( subcyc->SegListG0 );
	if ( subcyc->SegListG1 != (Segm_t *)NULL ) free( subcyc->SegListG1 );

	subcyc->SegListG0_siz = subcyc->totalAllSegs ;  subcyc->SegListG0_cnt = 0;
	subcyc->SegListG1_siz = subcyc->totalAllSegs ;  subcyc->SegListG1_cnt = 0;

	subcyc->SegListG0 = malloc(subcyc->SegListG0_siz * sizeof(Segm_t));
	subcyc->SegListG1 = malloc(subcyc->SegListG1_siz * sizeof(Segm_t));

	if ( subcyc->SegListG0==(Segm_t *)NULL || 
	     subcyc->SegListG1==(Segm_t *)NULL ) {
		puts ("Error (re)allocating memory in InitSegLists function");
		exit (1);
	}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

static void AddToSegSegList(Subcyc_t *subcyc, Segment_t *seg1, Segment_t *seg2, 
                            real8 dist2, int setSeg1Forces, int setSeg2Forces,
							int SubGroup)
{
        if (SubGroup == 0) {
		
			if (subcyc->SegSegListG0 == (SegSeg_t *)NULL || 
			    subcyc->SegSegListG0_cnt == subcyc->SegSegListG0_siz) {
				subcyc->SegSegListG0_siz = MAX(2*subcyc->SegSegListG0_siz, 10000);
				subcyc->SegSegListG0 = (SegSeg_t *)realloc(subcyc->SegSegListG0, 
				                        sizeof(SegSeg_t) * subcyc->SegSegListG0_siz);
				
				if ( subcyc->SegSegListG0 == (SegSeg_t *)NULL) {
					puts ("Error (re)allocating memory for SegSegListG0");
					exit (1);
				}
				
#ifdef _GPU_SUBCYCLE
				subcyc->SegSegListNodeInd = (int *)realloc(subcyc->SegSegListNodeInd, 
				                            sizeof(int) * 4 * subcyc->SegSegListG0_siz);
				subcyc->SegSegListArmInd  = (int *)realloc(subcyc->SegSegListArmInd, 
				                            sizeof(int) * 4 * subcyc->SegSegListG0_siz);
#endif
			}

			subcyc->SegSegListG0[subcyc->SegSegListG0_cnt].seg1          = seg1;
			subcyc->SegSegListG0[subcyc->SegSegListG0_cnt].seg2          = seg2;
			subcyc->SegSegListG0[subcyc->SegSegListG0_cnt].dist2         = dist2;
			subcyc->SegSegListG0[subcyc->SegSegListG0_cnt].flag          = 1;
			subcyc->SegSegListG0[subcyc->SegSegListG0_cnt].setSeg1Forces = setSeg1Forces;
			subcyc->SegSegListG0[subcyc->SegSegListG0_cnt].setSeg2Forces = setSeg2Forces;
			
#ifdef _GPU_SUBCYCLE
/*
 * 			These arrays are used to efficiently build the
 * 			reduction arrays when using the GPU subcycling.
 */
			subcyc->SegSegListNodeInd[4*subcyc->SegSegListG0_cnt+0] = seg1->node1->numInt++;
			subcyc->SegSegListNodeInd[4*subcyc->SegSegListG0_cnt+1] = seg1->node2->numInt++;
			subcyc->SegSegListNodeInd[4*subcyc->SegSegListG0_cnt+2] = seg2->node1->numInt++;
			subcyc->SegSegListNodeInd[4*subcyc->SegSegListG0_cnt+3] = seg2->node2->numInt++;
			
			subcyc->SegSegListArmInd[4*subcyc->SegSegListG0_cnt+0] = seg1->node1->armInt[seg1->armID12]++;
			subcyc->SegSegListArmInd[4*subcyc->SegSegListG0_cnt+1] = seg1->node2->armInt[seg1->armID21]++;
			subcyc->SegSegListArmInd[4*subcyc->SegSegListG0_cnt+2] = seg2->node1->armInt[seg2->armID12]++;
			subcyc->SegSegListArmInd[4*subcyc->SegSegListG0_cnt+3] = seg2->node2->armInt[seg2->armID21]++;
#endif

			subcyc->SegSegListG0_cnt += 1;
			
		} else if (SubGroup == 1) {
		
			if (subcyc->SegSegListG1 == (SegSeg_t *)NULL || 
			    subcyc->SegSegListG1_cnt == subcyc->SegSegListG1_siz) {
				subcyc->SegSegListG1_siz = MAX(2*subcyc->SegSegListG1_siz, 1000);
				subcyc->SegSegListG1 = (SegSeg_t *)realloc(subcyc->SegSegListG1, sizeof(SegSeg_t) *
				                                   subcyc->SegSegListG1_siz);
				
				if ( subcyc->SegSegListG1 == (SegSeg_t *)NULL) {
					puts ("Error (re)allocating memory for SegSegListG1");
					exit (1);
				}
			}

			subcyc->SegSegListG1[subcyc->SegSegListG1_cnt].seg1          = seg1;
			subcyc->SegSegListG1[subcyc->SegSegListG1_cnt].seg2          = seg2;
			subcyc->SegSegListG1[subcyc->SegSegListG1_cnt].dist2         = dist2;
			subcyc->SegSegListG1[subcyc->SegSegListG1_cnt].flag          = 1;
			subcyc->SegSegListG1[subcyc->SegSegListG1_cnt].setSeg1Forces = setSeg1Forces;
			subcyc->SegSegListG1[subcyc->SegSegListG1_cnt].setSeg2Forces = setSeg2Forces;

			subcyc->SegSegListG1_cnt += 1;
			
		} else if (SubGroup == 2) {
		
			if (subcyc->SegSegListG2 == (SegSeg_t *)NULL || 
			    subcyc->SegSegListG2_cnt == subcyc->SegSegListG2_siz) {
				subcyc->SegSegListG2_siz = MAX(2*subcyc->SegSegListG2_siz, 1000);
				subcyc->SegSegListG2 = (SegSeg_t *)realloc(subcyc->SegSegListG2, sizeof(SegSeg_t) *
				                                   subcyc->SegSegListG2_siz);
				
				if ( subcyc->SegSegListG2 == (SegSeg_t *)NULL) {
					puts ("Error (re)allocating memory for SegSegListG2");
					exit (1);
				}
			}

			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].seg1          = seg1;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].seg2          = seg2;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].dist2         = dist2;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].flag          = 1;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].setSeg1Forces = setSeg1Forces;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].setSeg2Forces = setSeg2Forces;

			subcyc->SegSegListG2_cnt += 1;
			
		} else if (SubGroup == 3) {
		
			if (subcyc->SegSegListG3 == (SegSeg_t *)NULL || 
			    subcyc->SegSegListG3_cnt == subcyc->SegSegListG3_siz) {
				subcyc->SegSegListG3_siz = MAX(2*subcyc->SegSegListG3_siz, 1000);
				subcyc->SegSegListG3 = (SegSeg_t *)realloc(subcyc->SegSegListG3, sizeof(SegSeg_t) *
				                                   subcyc->SegSegListG3_siz);
				
				if ( subcyc->SegSegListG3 == (SegSeg_t *)NULL) {
					puts ("Error (re)allocating memory for SegSegListG3");
					exit (1);
				}
			}

			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].seg1          = seg1;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].seg2          = seg2;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].dist2         = dist2;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].flag          = 1;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].setSeg1Forces = setSeg1Forces;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].setSeg2Forces = setSeg2Forces;

			subcyc->SegSegListG3_cnt += 1;
			
		} else if (SubGroup == 4) {
		
			if (subcyc->SegSegListG4 == (SegSeg_t *)NULL || 
			    subcyc->SegSegListG4_cnt == subcyc->SegSegListG4_siz) {
				subcyc->SegSegListG4_siz = MAX(2*subcyc->SegSegListG4_siz, 1000);
				subcyc->SegSegListG4 = (SegSeg_t *)realloc(subcyc->SegSegListG4, sizeof(SegSeg_t) *
				                                   subcyc->SegSegListG4_siz);
				
				if ( subcyc->SegSegListG4 == (SegSeg_t *)NULL) {
					puts ("Error (re)allocating memory for SegSegListG4");
					exit (1);
				}
			}

			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].seg1          = seg1;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].seg2          = seg2;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].dist2         = dist2;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].flag          = 1;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].setSeg1Forces = setSeg1Forces;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].setSeg2Forces = setSeg2Forces;

			subcyc->SegSegListG4_cnt += 1;
		
		}
		
        return;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void CellSegListMaker(Home_t *home)
{
		int         cellNum , cellID , armCount , cellSegCnt , i , arm;
		int         homeCells , homeNativeCells , homeDomain , numDomains;
		Subcyc_t   *subcyc;
		Param_t    *param;
		Cell_t     *cell;
        Node_t     *node, *nbr;
		
		homeCells       = home->cellCount;
		homeNativeCells = home->nativeCellCount;
		homeDomain      = home->myDomain;
		numDomains      = home->numDomains;
		param           = home->param;
		subcyc          = home->subcyc;

/*
 *		Allocate and initialize some temporary arrays we'll be needing.
 *
 *		For each cell native to or neighboring the current domain, an
 *		array of unique segments is built.  These arrays contain two classes
 *		of segments; "native" and "ghost" (see descriptions in the code
 *		below).  Any "native" segments in a cell's segment list will
 *		proceed "ghost" segments.  The lists are set up this way to
 *		simplify the process of insuring that forces on each pair
 *		of segments are only evaluated one time.
 *
 *		The cell segment arrays are set up as arrays of Segment_t
 *		structures,  each segment represented by a pair of node
 *		pointers and force component for each node of the segment.
 *		These force values will only represent the force on the
 *		nodes from the seg-seg force calcs done by this domain.
 *		These values will be summed with values from remote domains
 *		to get the final forces on the nodes/segments after all
 *		local calculations are done.
 */
		subcyc->totalNativeSegs = 0;
		
		if (subcyc->nativeSegCounts != (int *)NULL) free(subcyc->nativeSegCounts);
		if (subcyc->totalSegCounts  != (int *)NULL) free(subcyc->totalSegCounts );
		if (subcyc->cellSegLists    != (Segment_t **)NULL) {
			for (i = 0 ; i < subcyc->cellSegLists_siz ; i++) {
				if (subcyc->cellSegLists[i] == (Segment_t *)NULL) continue;
				free(subcyc->cellSegLists[i]);
			}
			free(subcyc->cellSegLists);
		}
		
		subcyc->nativeSegCounts = (int *)        calloc(1, sizeof(int)         * homeCells);
		subcyc->totalSegCounts  = (int *)        calloc(1, sizeof(int)         * homeCells);
		subcyc->cellSegLists    = (Segment_t **) calloc(1, sizeof(Segment_t *) * homeCells);

		for (cellNum = 0; cellNum < homeCells; cellNum++) {

			cellID = home->cellList[cellNum];
			cell   = home->cellKeys[cellID];

			if (cell->nodeCount == 0) continue;

/*
 *			Need to allocate a segment array large enough
 *			for all segments in the cell; could just do something
 *			like assume some maximum number of arms, multiply
 *			the node count by that factor and allocate space for
 *			that many pointer pairs... but for now do it the safe
 *			way and allocate 1 pointer pair per arm.
 */
			armCount = 0;
			node = cell->nodeQ;

			while (node != (Node_t *)NULL) {
				armCount += node->numNbrs;
				node = node->nextInCell;
			}
			
			subcyc->cellSegLists[cellNum] = (Segment_t *)calloc(1, sizeof(Segment_t) * armCount);
		}
		
		subcyc->cellSegLists_siz = homeCells;

/*
 *		Loop over all native cells adding "native" segments to the cell
 *		segment lists.  We only add native segments in this loop because
 *		all native segments in the array must proceed ghost segments (it
 *		makes things easier later on).  Ghost segments will be added
 *		to the arrays a little later.
 */
		for (cellNum = 0; cellNum < homeNativeCells; cellNum++) {

			cellID     = home->cellList[cellNum];
			cell       = home->cellKeys[cellID];
			node       = cell->nodeQ;
			cellSegCnt = 0;

/*
 *			Cell "native" segments are segments for which the dominant
 *			(i.e. owning) node of the segment is in the current cell and
 *			domain. 
 */
			for ( ; node != (Node_t *)NULL; node = node->nextInCell) {

				if (node->myTag.domainID != homeDomain) continue;

				for (arm = 0; arm < node->numNbrs; arm++) {

					nbr = GetNeighborNode(home, node, arm);

					if (nbr == (Node_t *)NULL) {
						printf("WARNING: Neighbor not found at %s line %d\n",
						        __FILE__, __LINE__);
						continue;
					}

					if (NodeOwnsSeg(home, node, nbr) == 0) continue;
					//if (OrderNodes(node, nbr) >= 0) continue;

					subcyc->cellSegLists[cellNum][cellSegCnt].node1    = node;
					subcyc->cellSegLists[cellNum][cellSegCnt].node2    = nbr;
					subcyc->cellSegLists[cellNum][cellSegCnt].armID12  = arm;
					subcyc->cellSegLists[cellNum][cellSegCnt].armID21  = GetArmID(home, nbr, node);
					subcyc->cellSegLists[cellNum][cellSegCnt].subGroup = 0;

					cellSegCnt++;
				}
			}

			subcyc->nativeSegCounts[cellNum] = cellSegCnt;
			subcyc->totalSegCounts [cellNum] = cellSegCnt;
			subcyc->totalNativeSegs         += cellSegCnt;
		}

/*
 *		Next add "ghost" segments to cell segment lists for
 *		all cells that are either native cells, or ghost cells
 *		which are NOT dominant over ALL cells native to the domain.
 *
 *		Note: A native cell may in fact only partially overlap
 *		the domain, hence the native cell may also contain ghost
 *		segments.
 *
 *		If there are NO native segments in this domain, there's
 *		no point in adding ghost segments to the lists because
 *		no force calcs will be done by this domain...
 */
		if (subcyc->totalNativeSegs == 0) {
			cellNum = homeCells;
		} else {
			cellNum = 0;
		}

		for (/* init cellNum above */ ; cellNum < homeCells; cellNum++) {

			cell       = home->cellKeys[home->cellList[cellNum]];
			node       = cell->nodeQ;
			cellSegCnt = subcyc->totalSegCounts[cellNum];

/*
 *			Cell "ghost" segments are comprised of segments owned by
 *			a "ghost" node.
 */
			for ( ; node != (Node_t *)NULL; node = node->nextInCell) {

				if (node->myTag.domainID == homeDomain) continue;

				for (arm = 0; arm < node->numNbrs; arm++) {

					nbr = GetNeighborNode(home, node, arm);

					if (nbr == (Node_t *)NULL) {
						printf("WARNING: Neighbor not found at %s line %d\n",
							    __FILE__, __LINE__);
						continue;
					}

					if (NodeOwnsSeg(home, node, nbr) == 0) continue;
				//	if (OrderNodes(node, nbr) >= 0) continue;

					subcyc->cellSegLists[cellNum][cellSegCnt].node1    = node;
					subcyc->cellSegLists[cellNum][cellSegCnt].node2    = nbr;
					subcyc->cellSegLists[cellNum][cellSegCnt].armID12  = arm;
					subcyc->cellSegLists[cellNum][cellSegCnt].armID21  = GetArmID(home, nbr, node);
					subcyc->cellSegLists[cellNum][cellSegCnt].subGroup = 0;
					cellSegCnt++;
				}
			}

			subcyc->totalSegCounts[cellNum] = cellSegCnt;
		}

/*
 *      Compute the total number of segments, ghost and native.
 */
        subcyc->totalAllSegs = 0;
        for (cellNum = 0; cellNum < homeCells; cellNum++) {
            subcyc->totalAllSegs += subcyc->totalSegCounts[cellNum];
        }

}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void SegSegListMaker(Home_t *home, int reqType)
{
		int         i , j , k , l , m;
		int         numNbrCells, cellNativeSegs , cellTotalSegs , SubAreaGroup;
		int         cellNum , cellID , armCount , cellSegCnt , nbrCellSegCnt , arm;
		int         homeCells , homeNativeCells , homeDomain , nbrCellID, numDomains;
		int         setSeg1Forces, setSeg2Forces, SubGroup   , SubThetaGroup, newNode;
		int         segGhost, seg1Ghost, seg2Ghost, nbrCellNatSegCnt;
		real8       dx , dy    , lgroup  , rgroup1S , rgroup3S , rgroupNN;
		real8       dz , dist2 , lgroupS , rgroup2S , rgroup4S;
		Param_t    *param;
		Subcyc_t   *subcyc;
		Cell_t     *cell, *nbrCell;
		Node_t     *node1 , *node2 , *node3 , *node4 , *nc1 , *nc2;
		Segment_t  *nbrSegList , *seg1 , *seg2;
		
		homeCells       = home->cellCount;
		homeNativeCells = home->nativeCellCount;
		homeDomain      = home->myDomain;
		numDomains      = home->numDomains;
		param           = home->param;
		subcyc          = home->subcyc;
		rgroup1S        = param->rg1 * param->rg1;
		rgroup2S        = param->rg2 * param->rg2;
		rgroup3S        = param->rg3 * param->rg3;
		rgroup4S        = param->rg4 * param->rg4;
		rgroupNN        = MAX(rgroup1S , rgroup2S) * (2.0 * 2.0);
		rgroupNN        = rgroup4S;
		lgroup          = param->minSeg * 0.25; 
		lgroupS         = lgroup * lgroup;
		
		TimerStart(home, SEG_LIST_MAKER);
		
/*
 *
 */
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			for (j = 0 ; j < 5 ; j++)
				node1->CommSend[j] = 0;
			
#ifdef _GPU_SUBCYCLE
			node1->numInt = 0;
			node1->armid = (int*)malloc(sizeof(int)*node1->numNbrs);
			node1->armInt = (int*)malloc(sizeof(int)*node1->numNbrs);
			for (j = 0; j < node1->numNbrs; j++) {
				node1->armInt[j] = 0;
			}
#endif
		}
		
/*
 *		Constructing the segment list corresponding to each cell of this domain
 */
		CellSegListMaker(home);

		InitSegLists(home);

/*
 *		Now go through and build a list of all the segment pairs
 *		for which forces need to be calculated.  (Note: in a segment
 *		pair, it's possible we don't need to get the forces for one
 *		of the segments in the pair.)
 *
 *		Note: The list of the native segments for which we have to calculate 
 *		forces will be built later,
 *
 *		Since the force calculation returns forces for all four nodes
 *		defining a pair of segments, we have to be sure we only
 *		compute the forces once for every pair of segments.  Additionally
 *		we don't need to compute forces between segment pairs if neither
 *		of the segments is native (i.e. has a node native to the current
 *		domain.)
 *
 *		The outer loop here only goes through the cells native to the
 *		domain because none of the other cells will have native segments.
 */
		
		subcyc->SegSegListG0_cnt = 0;
		subcyc->SegSegListG1_cnt = 0;
		subcyc->SegSegListG2_cnt = 0;
		subcyc->SegSegListG3_cnt = 0;
		subcyc->SegSegListG4_cnt = 0;
		
		if (param->elasticinteraction) {

		for (i = 0; i < homeNativeCells; i++) {
			
			cellNativeSegs = subcyc->nativeSegCounts[i];
			cellTotalSegs  = subcyc->totalSegCounts [i];
			
			cellID = home->cellList[i];
			cell   = home->cellKeys[cellID];

/*
 *			We'll need seg/seg forces between every native
 *			segment in this cell and all other segments following
 *			the segment in the cell's segment list. (Unless the
 *			the 2nd segment is a ghost segment, and the node owning
 *			the second segment is lower priority (force calc will
 *			then be done by domain owning segment 2)).
 */
			for (j = 0; j < cellNativeSegs; j++) {

				setSeg1Forces = 1;

				node1 = subcyc->cellSegLists[i][j].node1;
				node2 = subcyc->cellSegLists[i][j].node2;

/*
 *				If we're only doing a partial force calc, we don't
 *				need to update the forces for a native segment if 
 *				is not attached to a node marked for update.
 */
				if (reqType == PARTIAL) {
					if (((node1->flags & NODE_RESET_FORCES) == 0) &&
						((node2->flags & NODE_RESET_FORCES) == 0)) {
						setSeg1Forces = 0;
					}
				}

/*
 *				Now for segment pairs for which interactions must
 *				be computed.
 */
				for (k = j + 1; k < cellTotalSegs; k++) {

					setSeg2Forces = 1;

					node3 = subcyc->cellSegLists[i][k].node1;
					node4 = subcyc->cellSegLists[i][k].node2;

/*
 *					If we're only doing a partial force calc, we won't
 *					need forces for segment 2 if it is not attached to a node
 *					marked for update.
 */
					if (reqType == PARTIAL) {
						if (((node3->flags & NODE_RESET_FORCES) == 0) &&
						    ((node4->flags & NODE_RESET_FORCES) == 0)) {
							setSeg2Forces = 0;
						}
					}

/*
 *					If neither of the segments needs forces updated, skip it.
 */
					if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) continue;

/*
 *					If segment 2 is a ghost, only do the force calc if the
 *					node owning segment 1 has lower priority than the
 *					node owning segment 2. (this if statement becomes true only
 *					when some part of the cell is in the current domain and 
 *					some parts are out of domain).
 *					segGhost = 1 : this segment pair will be done by another domain
 *					segGhost = 0 : this segment pair will be done by this domain
*/
					if ((k >= cellNativeSegs) && (NodeOwnsSeg(home, node1, node3))) 
						segGhost = 1;
					else
						segGhost = 0;
					
/*
 *					This segment pair should be placed in either group 0 or 1 
 *					depending on the distance between them. So first calculate 
 *					the distance. In the mean time check if the line-tension 
 *					force calculation of the segments should be placed in group 1.
 */
					GetMinDistSegSeg( home, node1, node2, node3, node4, &dist2);
					if (dist2 < 0) continue;
/*
 *					Assign each seg-seg interaction to subgroup based on their distance
 */
					newNode = 0;
					if ((node1->newNode > 0) || (node2->newNode > 0) ||
						(node3->newNode > 0) || (node4->newNode > 0) ) newNode = 1;
						
					if (dist2 < rgroup1S) {
						SubGroup = 1;
					} else if (dist2 < rgroup2S ) {
						SubGroup = 2;
					} else if (dist2 < rgroup3S ) {
						SubGroup = 3;
					} else if (dist2 < rgroup4S || (newNode && dist2 < rgroupNN)) {
						SubGroup = 4;
					} else {
						SubGroup = 0;
					}
/*
 *					This segment pair should be placed in one of the groups 0-4 
 *					Add segment pair to segment pair list, unless the force
 *					will be computed by another domain in which case
 *					we add the segment to the ghost list.
 */
					if (param->sendSubGroupForc==0 && SubGroup > 0) continue;
					
                    if (segGhost) {
						node1->CommSend[SubGroup] = 1;
						node2->CommSend[SubGroup] = 1;
						node3->CommSend[SubGroup] = 1;
						node4->CommSend[SubGroup] = 1;
                    } else {
						
						AddToSegSegList(subcyc, &(subcyc->cellSegLists[i][j]), 
							                    &(subcyc->cellSegLists[i][k]), 
						                        dist2 , setSeg1Forces, 
												setSeg2Forces, SubGroup);
                    }
				}

			}  /* Loop over native segments */
			
/*
 *			Next loop over all the neighbours of the current native cell. If the 
 *			current cell has priority over the the neighbouring cell we do need 
 *			force calcs between these pairs in this loop; the segment pair will 
 *			either be handled by one of the other iterations of this loop, or 
 *			by the remote domain owning the segments in the other cell.
 *
 *			Note: Cell ids used here convert to ranges from zero -> num[XYZ]cells+1 
 *			allowing for periodic cells. But nbrCellID gets converted to the 
 *			base cell index if it is a periodic cell.
 */

			numNbrCells = cell->nbrCount;

			for (j = 0; j < numNbrCells; j++) {

				nbrCellID = cell->nbrList[j];
				nbrCell = home->cellKeys[nbrCellID];

				if (nbrCell->baseIdx >= 0) nbrCellID = nbrCell->baseIdx;

/*
 *              If the neighbour cell has priority over the
 *              current cell we need to calculate seg/seg forces
 *              between native segs in the current cell and the
 *              segments in the neighbouring cell.
 */
				if (CellPriority(home, cellID, nbrCellID) >= 0)
					seg1Ghost = 1;
				else
					seg1Ghost = 0;
					

                for (k = 0; k < homeCells; k++) {
                    if (nbrCellID == home->cellList[k]) break;
                }

                nbrSegList       = subcyc->cellSegLists[k];
                nbrCellSegCnt    = subcyc->totalSegCounts[k];
				nbrCellNatSegCnt = subcyc->nativeSegCounts[k];

/*
 *              If there are no segments in the neighbouring cell, no
 *              need to do anything more with this neighbour cell.
 */
                if (nbrCellSegCnt < 1) continue;

                for (k = 0; k < cellNativeSegs; k++) {

                    node1 = subcyc->cellSegLists[i][k].node1;
                    node2 = subcyc->cellSegLists[i][k].node2;

                    setSeg1Forces = 1;

/*
 *                  If we're only doing a partial force calc, we don't
 *                  need forces for segment 1 if it is not attached to a node
 *                  marked for update.
 */
                    if (reqType == PARTIAL) {
                        if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                            ((node2->flags & NODE_RESET_FORCES) == 0)) {
                            setSeg1Forces = 0;
                        }
                    }

                    for (l = 0; l < nbrCellSegCnt; l++) {

                        node3 = nbrSegList[l].node1;
                        node4 = nbrSegList[l].node2;
						
						if (seg1Ghost == 1) {
							if (node3->myTag.domainID == homeDomain &&
							    node4->myTag.domainID == homeDomain ) continue;
						}

                        setSeg2Forces = 1;

/*
 *                      If we're only doing a partial force calc, we don't
 *                      need forces for segment 2 if it is not attached to
 *                      a node marked for update.
 */
						if (reqType == PARTIAL) {
							if (((node3->flags & NODE_RESET_FORCES) == 0) &&
							    ((node4->flags & NODE_RESET_FORCES) == 0)) {
								setSeg2Forces = 0;
							}
						}

/*
 *						If the 2nd segment is native, we probably need the 
 *						forces, but if segment 2 is a ghost, only do the
 *						calc if the node owning segment 1 has lower priority
 *						than the node owning segment 2.
 */
                        if ((l >= nbrCellNatSegCnt) && (NodeOwnsSeg(home, node1, node3)))
							seg2Ghost = 1;
						else
							seg2Ghost = 0;

						if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) continue;

/*
 *						This segment pair should be placed in one of the groups 0-4 
 *						depending on the distance between them. So first calculate 
 *						the distance.
 */
						GetMinDistSegSeg( home, node1, node2, node3, node4, &dist2);
						if (dist2 < 0) continue;
/*
 *						Assign each seg-seg interaction to subgroup based on their distance
 */
						newNode = 0;
						if ((node1->newNode > 0) || (node2->newNode > 0) ||
							(node3->newNode > 0) || (node4->newNode > 0) ) newNode = 1;
						
						if (dist2 < rgroup1S) {
							SubGroup = 1;
						} else if (dist2 < rgroup2S ) {
							SubGroup = 2;
						} else if (dist2 < rgroup3S ) {
							SubGroup = 3;
						} else if (dist2 < rgroup4S || (newNode && dist2 < rgroupNN)) {
							SubGroup = 4;
						} else {
							SubGroup = 0;
						}
						//	SubGroup = 1;
/*
 *                      Add segment pair to segment pair list, unless the force
 *                      will be computed by another domain in which case
 *                      we add the segment to the ghost list.
 */
						if (param->sendSubGroupForc==0 && SubGroup > 0) continue;
						
						if (seg1Ghost || seg2Ghost) {
							node1->CommSend[SubGroup] = 1;
							node2->CommSend[SubGroup] = 1;
							node3->CommSend[SubGroup] = 1;
							node4->CommSend[SubGroup] = 1;
                        } else {
							
							AddToSegSegList(subcyc, &(subcyc->cellSegLists[i][k]), &nbrSegList[l], 
											dist2 , setSeg1Forces, setSeg2Forces, SubGroup);
                        }
					}
				}
			}  /* for (j = 0; j < numNbrCells...) */
			
		} /* for (i = 0; i < homeCells...) */
		
		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		if (param->sendSubGroupForc == 0) {
			int lStart;
			
			for (i = 0; i < homeCells; i++) {
				
				cellNativeSegs = subcyc->nativeSegCounts[i];
				cellTotalSegs  = subcyc->totalSegCounts [i];
				
				for (j = i; j < homeCells; j++) {
					nbrCellNatSegCnt = subcyc->nativeSegCounts[j];
					nbrCellSegCnt    = subcyc->totalSegCounts [j];
					
					for (k = 0; k < cellTotalSegs; k++) {

						setSeg1Forces = 1;

						node1 = subcyc->cellSegLists[i][k].node1;
						node2 = subcyc->cellSegLists[i][k].node2;

/*
 *						If we're only doing a partial force calc, we don't
 *						need to update the forces for a native segment if 
 *						is not attached to a node marked for update.
 */
						if (reqType == PARTIAL) {
							if (((node1->flags & NODE_RESET_FORCES) == 0) &&
								((node2->flags & NODE_RESET_FORCES) == 0)) {
								setSeg1Forces = 0;
							}
						}
						
						if (i == j) lStart = k+1;
						else        lStart = 0  ;
						
						for (l = lStart; l < nbrCellSegCnt; l++) {

							setSeg2Forces = 1;

							node3 = subcyc->cellSegLists[j][l].node1;
							node4 = subcyc->cellSegLists[j][l].node2;

							if (reqType == PARTIAL) {
								if (((node3->flags & NODE_RESET_FORCES) == 0) &&
									((node4->flags & NODE_RESET_FORCES) == 0)) {
									setSeg2Forces = 0;
								}
							}
							
/*
 *							If neither of the segments needs forces updated, skip it.
 */
							if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) continue;
							
/*
 *							If none of the nodes are in this domain, skip it
 */
							if (node1->myTag.domainID != homeDomain &&
								node2->myTag.domainID != homeDomain &&
								node3->myTag.domainID != homeDomain &&
								node4->myTag.domainID != homeDomain ) continue;
								
/*
 *							If one of the nodes is ghost, mark all of the four nodes
 *							for velcoity/position communication.
 */
							segGhost = 0;
							if (node1->myTag.domainID != homeDomain ||
								node2->myTag.domainID != homeDomain ||
								node3->myTag.domainID != homeDomain ||
								node4->myTag.domainID != homeDomain ) segGhost = 1;
/*
 *							This segment pair should be placed in either group 0-4 
 *							depending on the distance between them. So first calculate 
 *							the distance.
 */
							GetMinDistSegSeg( home, node1, node2, node3, node4, &dist2);
							if (dist2 < 0) continue;
/*
 *							Assign each seg-seg interaction to subgroup based on their distance
 */
							newNode = 0;
							if ((node1->newNode > 0) || (node2->newNode > 0) ||
								(node3->newNode > 0) || (node4->newNode > 0) ) newNode = 1;
								
							if (dist2 < rgroup1S) {
								SubGroup = 1;
							} else if (dist2 < rgroup2S ) {
								SubGroup = 2;
							} else if (dist2 < rgroup3S ) {
								SubGroup = 3;
							} else if (dist2 < rgroup4S || (newNode && dist2 < rgroupNN)) {
								SubGroup = 4;
							} else {
								SubGroup = 0;
							}
							
							if (SubGroup == 0) continue;
/*
 *                 	    	Add segment pair to segment pair list
 */
							AddToSegSegList(subcyc, &(subcyc->cellSegLists[i][k]), 
													&(subcyc->cellSegLists[j][l]), 
													dist2 , setSeg1Forces, 
													setSeg2Forces, SubGroup);
							
							if (segGhost) {
								node1->CommSend[SubGroup] = 1;
								node2->CommSend[SubGroup] = 1;
								node3->CommSend[SubGroup] = 1;
								node4->CommSend[SubGroup] = 1;
							}
						}
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////
		} /* param->elasticinteraction */
		
/*
 *		We have built the list of segments pair. Now go and generate the 
 *		segment lists. These lists are used for the calculations of self-forces
 *		and forces due to external stresses, etc.
 */		
		for (i = 0; i < homeCells; i++) {
			if (param->sendSubGroupForc && i >= homeNativeCells) continue;

			for (j = 0; j < subcyc->totalSegCounts[i]; j++) {
				if (param->sendSubGroupForc && j >= subcyc->nativeSegCounts[i]) continue;

				setSeg1Forces = 1;

				node1 = subcyc->cellSegLists[i][j].node1;
				node2 = subcyc->cellSegLists[i][j].node2;

/*
 *              If neither node is owned by this domain, skip this segment.
 */
                if ((node1->myTag.domainID != homeDomain) && 
                    (node2->myTag.domainID != homeDomain)) continue;

/*
 *              If one node is not owned by this domain, flag the segment
 *              as a ghost.
 */
				segGhost = 0;
				if ((node1->myTag.domainID != homeDomain) || 
				    (node2->myTag.domainID != homeDomain)) segGhost = 1;

/*
 *				If we're only doing a partial force calc, we don't
 *				need to update the forces for a native segment if 
 *				is not attached to a node marked for update.
 */
				if (reqType == PARTIAL) {
					if (((node1->flags & NODE_RESET_FORCES) == 0) &&
						((node2->flags & NODE_RESET_FORCES) == 0)) {
						setSeg1Forces = 0;
					}
				}

/*
 *				If necessary, add native segment to list of G0 and G1 segs
 *				for which we need to calculate forces from interactions
 *				other than seg/seg interactions (i.e. self force,
 *				osmotic force, remote force, etc.)
 */
				if (setSeg1Forces) {
					subcyc->SegListG1[subcyc->SegListG1_cnt].seg  = &subcyc->cellSegLists[i][j];
					subcyc->SegListG1[subcyc->SegListG1_cnt].flag = 1;
					subcyc->SegListG1_cnt++;
				}
                
/*
 *				If this is a ghost, add it to the ghost segment list so that
 *				it gets flagged for parallel communication.
 */
				if (segGhost) {
					node1->CommSend[1] = 1;
					node2->CommSend[1] = 1;
				}

			}
		}

		
/*
 *		The interactions within the enlarged cutoff of the 'new nodes' were 
 *		placed in G4. So turn of all the newNode flags.
 */
		for (i = 0 ; i < home->newNodeKeyPtr ; i++) {
            if ( (node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			node1->newNode = 0;
		}
		
		TimerStop(home, SEG_LIST_MAKER);
}
