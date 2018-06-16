#ifdef _GPU_SUBCYCLE

#if PARALLEL
#error Cannot compile the GPU version in parallel
#endif

#ifndef _SubCycGPU_h
#define _SubCycGPU_h

typedef struct {
	int  deviceID;
	int  maxBlocks;
	int  maxThreadsPerBlock;
} DeviceProp_t;

typedef struct {
	int     nArms, nodeCount, segPairCnt;
	int    *n12, *n34;
	double *cc, *b1, *b2;
	double *fpair;
} SplitSegSeg_t;

#ifdef __cplusplus

extern "C" {
	void RemoteSigbSub(Home_t *home);
	void InitializeParadisGPU(Home_t *home);
	void SubcycleIntegratorGPU(Home_t *home);
	void LocalSegForcesGPU(Home_t *home, int segPairListCnt, SegmentPair_t *segPairList, double *fpair);
	void SetOneNodeForceGPU(Home_t *home, SplitSegSeg_t *splitSegSegList);
	void InitializeNodeForceGPU(Home_t *home);
	void FinalizeNodeForceGPU(Home_t *home);
#ifdef __NVCC__
	void SelectCudaDevice(int deviceID, DeviceProp_t *deviceProp);
	int GetCudaCores(cudaDeviceProp devProp);
	void GetThreadsPerBlock(DeviceProp_t *deviceProp, int maxThreadsSize, int *threadsPerBlock, int *blockSize);
#endif
}

#ifdef __NVCC__
typedef struct {
	
	int      threadsPerBlock;
	
	int      nodeCount;
	int      segCount;
	int      armCount;
	int      blocksNodes;
	int      blocksSegs;
	int      blocksArms;
	
	double   *esig;
	int      *n;
	double3  *r, *r0, *b, *f, *v, *v0, *rkf;
	double   *mob, *B;
	double3  *fmm, *cc;
	
	double   *e1, *e2;
	thrust::device_ptr<double> e1_ptr;
	thrust::device_ptr<double> e2_ptr;
	
	int2     *g1pos;
	int      *g1ind;
	int      *g1arms;
	int2     *g1;
	double3  *fseg, *farms;
	double   *Bseg;
	
	int      nSegSeg0, blocksSegSegs0;
	int2     *g0pos;
	int      *g0ind;
	int2     *g0arms_pos;
	int      *g0arms_ind;
	int2     *g0;
	double3  *f0;
	int      *g0flag;
	double   *g0dist2;
	
	int      nSegSeg2, blocksSegSegs2;
	int2     *g2pos;
	int      *g2ind;
	int2     *g2arms_pos;
	int      *g2arms_ind;
	int2     *g2;
	double3  *f2;
	
	int      nSegSeg3, blocksSegSegs3;
	int2     *g3pos;
	int      *g3ind;
	int2     *g3arms_pos;
	int      *g3arms_ind;
	int2     *g3;
	double3  *f3;
	
	int      nSegSeg4, blocksSegSegs4;
	int2     *g4pos;
	int      *g4ind;
	int2     *g4arms_pos;
	int      *g4arms_ind;
	int2     *g4;
	double3  *f4;
	
	int4     *pair;
	double   *fpair;
	double3  *b1, *b2;
	
} Device_t;
#endif

#else
void RemoteSigbSub(Home_t *home);
void InitializeParadisGPU(Home_t *home);
void SubcycleIntegratorGPU(Home_t *home);
void LocalSegForcesGPU(Home_t *home, int segPairListCnt, SegmentPair_t *segPairList, double *fpair);
void SetOneNodeForceGPU(Home_t *home, SplitSegSeg_t *splitSegSegList);
void InitializeNodeForceGPU(Home_t *home);
void FinalizeNodeForceGPU(Home_t *home);
#endif

#endif
#endif
