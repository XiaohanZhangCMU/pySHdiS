#ifndef __SH_H__
#define __SH_H__

#ifdef _SHDIS

int GetNumPhysicNodes(Home_t *home);
#if 0
void GetNodeList(Home_t *home, double* XXX, double* YYY, double* ZZZ);
#else
void GetNodeList(Home_t *home, double* XYZ);
#endif
void SH_calc_stress(Home_t *home, double* stress);
void ComputeSH1SegSigbRem(Home_t *home,
			  Node_t *node, Node_t *nbr,
			  int armID1, int armID2);

#endif
#endif // __SH_H__
