/*--------------------------------------------------------------------------
 *
 *	Node.h	Define the struct that holds all relevant data for a single
 *		node, either native or ghost
 *
 *		Notice: make depend if new files are added, otherwise always
 *			make clean after changing this file  Wei Cai 04/09/2002
 *
 *------------------------------------------------------------------------*/

#ifndef _SubCyc_h
#define _SubCyc_h


typedef struct {
	int        forcesSet, subGroup;
	int        armID12  , armID21;
	real8      f1[3]    , f2[3]  ;
	Node_t     *node1   , *node2 ;
#ifdef _OPENMP
	omp_lock_t segLock;
#endif
#ifdef _GPU_SUBCYCLE
	int     subindex;
#endif
} Segment_t;

typedef struct {
	Segment_t   *seg;
	int          flag;
} Segm_t;

typedef struct {
	Segment_t   *seg1;
	Segment_t   *seg2;
	int          setSeg1Forces;
	int          setSeg2Forces;
	int          flag;
	real8        dist2;
} SegSeg_t;

typedef struct {
	Node_t    ***ShortRange_SegSegList;
	int          Size_SegSegList , max_SegSegList;
	
	int          numSubCycle1, numSubCycle2, numSubCycle3, numSubCycle4;
	real8        Group1Frac  , Group2Frac  , Group3Frac  , Group4Frac  ;
	
	int          SegSegListG0_siz , SegSegListG0_cnt , SegListG0_siz , SegListG0_cnt;
	int          SegSegListG1_siz , SegSegListG1_cnt , SegListG1_siz , SegListG1_cnt;
	int          SegSegListG2_siz , SegSegListG2_cnt ;
	int          SegSegListG3_siz , SegSegListG3_cnt ;
	int          SegSegListG4_siz , SegSegListG4_cnt ;
	
	Segm_t      *SegListG0        , *SegListG1       ;
	SegSeg_t    *SegSegListG0     , *SegSegListG1    , *SegSegListG2 ,  
				*SegSegListG3     , *SegSegListG4    ;
				
	Segment_t  **cellSegLists;
	int         *totalSegCounts   , *nativeSegCounts ;
	int          totalNativeSegs  ,  totalAllSegs    ;
	int          cellSegLists_siz ;
	
	real8      **sigbFMM;
	
#ifdef _GPU_SUBCYCLE
	int         *SegSegListNodeInd;
	int         *SegSegListArmInd;
#endif
	
} Subcyc_t;

//typedef struct _subcyc ;

#endif
