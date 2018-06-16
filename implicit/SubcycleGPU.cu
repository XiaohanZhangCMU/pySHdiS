/*-------------------------------------------------------------------------
 *
 *	This file contains the GPU implementation of the subcycling
 * 	time-integrator algorithm.
 * 
 * 	As of now, the GPU implementation has the following limitations:
 * 		- 	it can only be run with 1 CPU (serial mode)
 * 		- 	it can only be used with the FMM and cannot be used 
 * 			with the Rijm table
 * 		-	it can only be used with the FCC_0 mobility law
 * 		- 	it can not be used with rotated frames (must set useLabFrame = 0)
 * 		-	it cannot be used with interactions in subgroup 1
 * 			(i.e. one must set rg1 = 0)
 * 
 * 	Nicolas Bertin, 06/27/2017
 *
 *-----------------------------------------------------------------------*/

#ifdef _GPU_SUBCYCLE

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <cuda.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>

#include "Util.h"
#include "Home.h"
#include "SubcycleGPU.h"
#include "Mobility.h"

/*------------------------------------------------------------------------
 *
 *      Function:    HandleErrorGPU
 *
 *-----------------------------------------------------------------------*/
#define HANDLE_ERROR(err) (HandleErrorGPU(err, __FILE__, __LINE__ ))
static void HandleErrorGPU(cudaError_t err, const char *file, int line) {
	if (err != cudaSuccess) {
		printf( "%s in %s at line %d\n", cudaGetErrorString(err), file, line);
		exit(EXIT_FAILURE);
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    CheckErrorGPU
 *
 *-----------------------------------------------------------------------*/
void CheckErrorGPU(const char *message) {
	cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess) {
		fprintf(stderr,"ERROR %s: %s\n", message, cudaGetErrorString(error));
		exit(-1);
	}                        
}

/*------------------------------------------------------------------------
 *
 *      Function:    SetVariablesGPU
 *
 *-----------------------------------------------------------------------*/
__device__ double  MU, NU, a;
__device__ double  TensionFactor;
__device__ double  MobEdge, MobScrew;
__device__ double  Lx, Ly, Lz;
__device__ double  invLx, invLy, invLz;
__device__ double  rTol, rTolth, rTolrel;
__device__ double  coreEnergy;
__device__ double3 boxc;

__global__ void SetVariablesGPU(double _MU, double _NU, double _a, double _tf, 
                                double _MobEdge, double _MobScrew,
                                double _Lx, double _Ly, double _Lz,
                                double _xc, double _yc, double _zc,
                                double _rTol, double _rTolth, double _rTolrel, double _Ecore)
{
	MU = _MU;
	NU = _NU;
	a = _a;
	TensionFactor = _tf;
	MobEdge = _MobEdge;
	MobScrew = _MobScrew;
	Lx = _Lx;
	Ly = _Ly;
	Lz = _Lz;
	if (_Lx == 0.0) invLx = 0.0;
	else invLx = 1.0 / _Lx;
	if (_Ly == 0.0) invLy = 0.0;
	else invLy = 1.0 / _Ly;
	if (_Lz == 0.0) invLz = 0.0;
	else invLz = 1.0 / _Lz;
	boxc.x = _xc;
	boxc.y = _yc;
	boxc.z = _zc;
	rTol = _rTol;
	rTolth = _rTolth;
	rTolrel = _rTolrel;
	coreEnergy = _Ecore;
}

__device__ void SegSegForceIsotropicGPU(double3 r1, double3 r2, double3 r3, double3 r4, double3 b1, double3 b2,
                                        double3 &fn1, double3 &fn2, double3 &fn3, double3 &fn4);
__device__ void SegSegForceIsotropicCorrGPU(double3 r1, double3 r2, double3 r3, double3 r4, double3 b1, double3 b2,
                                            double3 &fn1, double3 &fn2, double3 &fn3, double3 &fn4);
__device__ void SpecialSegSegForceHalfGPU(double3 r1, double3 r2, double3 r3, double3 r4, 
                                          double3 b1, double3 b2, double3 &fn3, double3 &fn4);

/*------------------------------------------------------------------------
 *
 *      Function:    SpecialSegSegForceHalfGPU
 *
 *-----------------------------------------------------------------------*/
__device__ void SpecialSegSegForceHalfGPU(double3 r1, double3 r2, double3 r3, double3 r4, 
                                          double3 b1, double3 b2, double3 &fn3, double3 &fn4)
{
        int i, j , alt1[3]={1,2,0}, alt2[3]={2,0,1};
        real8 eps, ecrit, c, a2, d2, a2_d2, a2d2inv;
        real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
        real8 f3[3], f4[3];
        real8 vec1[3], vec2[3], t[3], nd[3];
        real8 temp1;
        real8 R[3], Rdt, x1mod[3], x2mod[3];
        real8 oneoverL;
        real8 y[2], z[2], yv[4], zv[4], ypz[4], ymz[4];
        real8 Ra[4], Rainv[4], Log_Ra_ypz[4];
        real8 temp, tmp[8];
        real8 common1[4], common2[3], common3[3];
        real8 magdiff, diffMag2, x1modMag2, x2modMag2;
        real8 f_003v[4], f_103v[4], f_113v[4], f_213v[4];
        real8 f_005v[4], f_105v[4], f_115v[4], f_215v[4];
        real8 f_003, f_103, f_113, f_213;
        real8 f_005, f_105, f_115, f_215;
        real8 Fint_003, Fint_113, Fint_005, Fint_115;
        real8 I_003[3], I_113[3], I_005[3], I_115[3];
        real8 m4p, m8p, m4pn, a2m4pn, a2m8p;
        real8 tdb, tdbp, nddb, bpctdb, bpctdnd;
        real8 bct[3], bpct[3], ndct[3], bpctct[3];
        real8 cotanthetac;
        real8 pivalue=3.141592653589793;

		ecrit  = 1e-4;
        cotanthetac = sqrt((1 - ecrit*1.01) / (ecrit*1.01));
        
        eps    = 1e-12;
        a2     = a*a;
        m4p    = 0.25 * MU / pivalue;
        m8p    = 0.5 * m4p;
        m4pn   = m4p / ( 1 - NU );
        a2m4pn = a2 * m4pn;
        a2m8p  = a2 * m8p;
            
        fn3.x = 0.0;
        fn3.y = 0.0;
        fn3.z = 0.0;
         
        fn4.x = 0.0;
        fn4.y = 0.0;
        fn4.z = 0.0;
        
        x1[0]=r1.x;
        x1[1]=r1.y;
        x1[2]=r1.z;
        x2[0]=r2.x;
        x2[1]=r2.y;
        x2[2]=r2.z;
        x3[0]=r3.x;
        x3[1]=r3.y;
        x3[2]=r3.z;
        x4[0]=r4.x;
        x4[1]=r4.y;
        x4[2]=r4.z;
        
        b[0]=b2.x;
        b[1]=b2.y;
        b[2]=b2.z;
        bp[0]=b1.x;
        bp[1]=b1.y;
        bp[2]=b1.z;
        
        #pragma unroll
        for(i=0;i<3;i++) { 
            vec1[i]=x4[i]-x3[i];
            vec2[i]=x2[i]-x1[i];
        }

        temp1=0.0e0;
        
        #pragma unroll
        for(i=0;i<3;i++) { 
            temp1+=vec1[i]*vec1[i];
        }

        oneoverL =1/sqrt(temp1);
        
        for(i=0;i<3;i++) { 
            t[i]=vec1[i]*oneoverL;
        }
        
        c=0.0e0;
        
        #pragma unroll
        for(i=0;i<3;i++) { 
            c+=t[i]*vec2[i];
        }

        if (c < 0) {
			#pragma unroll
            for(i=0;i<3;i++) { 
                temp=x2[i];
                x2[i]=x1[i];
                x1[i]=temp;
                bp[i]=-bp[i];
                vec2[i]=-vec2[i];
            }         
        }
/*
 *      Find f3 and f4, but only if at least one of the segment
 *      endpoints is local to the domain.
 */
        temp=0.0e0;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            temp+=vec2[i]*t[i];
        }
        
        #pragma unroll
        for (i=0;i<3;i++) {
            x2mod[i]=x1[i]+temp*t[i];
        }
        
        #pragma unroll
        for (i=0;i<3;i++) {
            vec2[i]=x2[i]-x2mod[i];
        }
               
        temp=0.0e0;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            temp+=vec2[i]*vec2[i];
        }
            
        magdiff=sqrt(temp);
        temp=magdiff*0.5e0 * cotanthetac;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            vec1[i]=temp*t[i];
        }
        
        #pragma unroll
        for (i=0;i<3;i++) {
            x1mod[i]=x1[i]+0.5e0*vec2[i]+vec1[i];
            x2mod[i]+=0.5e0*vec2[i]-vec1[i];
        }
        
        #pragma unroll
        for (i=0;i<3;i++) {
            R[i]=0.5e0*((x3[i]+x4[i])-(x1mod[i]+x2mod[i]));
        }
        
        Rdt=0.0e0;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            Rdt+=R[i]*t[i];
        }
        
        #pragma unroll
        for (i=0;i<3;i++) {
            nd[i]=R[i]-Rdt*t[i];
        }
        
        d2=0.0e0;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            d2+=nd[i]*nd[i];
        }    
        
        #pragma unroll
        for (j=0;j<2;j++) {
            y[j]=0.0e0;
            z[j]=0.0e0;
        }  
        
        #pragma unroll
        for (i=0;i<3;i++) {
            y[0]+=x3[i]*t[i];
            y[1]+=x4[i]*t[i];
            z[0]+=-x1mod[i]*t[i];
            z[1]+=-x2mod[i]*t[i];
        } 
        
        #pragma unroll
        for (j=0;j<2;j++) {
            yv[2*j]=y[j];
            yv[2*j+1]=y[j];
            zv[j]=z[j];
            zv[j+2]=z[j];
        }    
            
        a2_d2 = a2 + d2;   
        
        #pragma unroll
        for (j=0;j<4;j++) {
            ypz[j] = yv[j] + zv[j];
            ymz[j] = yv[j] - zv[j];
        }
        
        #pragma unroll
        for (j=0;j<4;j++) {
            tmp[j]=a2_d2 + ypz[j]*ypz[j];
        }
        
        #pragma unroll
        for (j=0;j<4;j++) {
            Ra[j]=sqrt(tmp[j]);
        }
        
        #pragma unroll
        for (j=0;j<4;j++) {
            Rainv[j]=1.0e0/Ra[j];
        }

        a2d2inv = 1.0e0 / a2_d2;
        
        #pragma unroll
        for (j=0;j<4;j++) {
            tmp[j]=Ra[j] + ypz[j];
			tmp[j+4]=Ra[j]-ypz[j];
        }
        
        #pragma unroll
        for (j=0;j<4;j++) {
            Log_Ra_ypz[j]=0.5e0*(log(tmp[j])-log(tmp[j+4]));
        }
        
        #pragma unroll
        for (j=0;j<4;j++) {
            common1[j] = ymz[j] * Ra[j] * a2d2inv;
            f_115v[j] = -a2d2inv * ypz[j] * Rainv[j];
        }
        
        temp=2.0e0*a2d2inv;
        
        #pragma unroll
        for (j=0;j<4;j++) {
            f_003v[j] = Ra[j];
            f_103v[j] = Log_Ra_ypz[j] - common1[j];
            f_113v[j] = -Log_Ra_ypz[j];
            f_213v[j] = zv[j]*Log_Ra_ypz[j] - Ra[j];
            f_005v[j] = temp*Ra[j] - Rainv[j];
            f_105v[j] = common1[j] - yv[j]*Rainv[j];
            f_215v[j] =  Rainv[j] - zv[j] * f_115v[j];
        }
        
        f_003 = 0.0e0;
        f_103 = 0.0e0;
        f_113 = 0.0e0;
        f_213 = 0.0e0;
        f_005 = 0.0e0;
        f_105 = 0.0e0;
        f_115 = 0.0e0;
        f_215 = 0.0e0;
        
        #pragma unroll
        for (j=1;j<3;j++) {
            f_003v[j] = -f_003v[j];
            f_103v[j] = -f_103v[j];
            f_113v[j] = -f_113v[j];
            f_213v[j] = -f_213v[j];
            f_005v[j] = -f_005v[j];
            f_105v[j] = -f_105v[j];
            f_115v[j] = -f_115v[j];
            f_215v[j] = -f_215v[j];
        }
        
        #pragma unroll
        for (j=0;j<4;j++) {
            f_003 += f_003v[j];
            f_103 += f_103v[j];
            f_113 += f_113v[j];
            f_213 += f_213v[j];
            f_005 += f_005v[j];
            f_105 += f_105v[j];
            f_115 += f_115v[j];
            f_215 += f_215v[j];
        }

        f_103 *= -0.5e0;    
        f_003 *=  a2d2inv;
        f_005 *=  a2d2inv;
        f_105 *=  a2d2inv;  
          
        #pragma unroll
        for (i=0;i<3;i++) {
            bct[i]=b[alt1[i]]*t[alt2[i]] - b[alt2[i]]*t[alt1[i]];
            bpct[i]=bp[alt1[i]]*t[alt2[i]] - bp[alt2[i]]*t[alt1[i]];
            ndct[i]=nd[alt1[i]]*t[alt2[i]] - nd[alt2[i]]*t[alt1[i]];
        }
        
        tdb=0.0e0;
        tdbp=0.0e0;
        nddb=0.0e0;
        bpctdb=0.0e0;
        bpctdnd=0.0e0;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            tdb += t[i]*b[i];
            tdbp+= t[i]*bp[i];
            nddb+= nd[i]*b[i];
            bpctdb += bpct[i]*b[i];
            bpctdnd += bpct[i]*nd[i];
            
        }
            
        temp = tdb*tdbp;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            bpctct[i] = tdbp*t[i] - bp[i];
            common2[i] = temp*nd[i];
            common3[i] = bpctdnd*bct[i];
        }   

        tmp[0]=(m4pn-m4p)*tdb;
        tmp[1]=m4pn*bpctdnd*nddb;
        tmp[2]=a2m8p*tdb;
        tmp[3]=m4pn*bpctdnd*tdb;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            I_003[i] = m4pn*(nddb*bpctct[i] + bpctdb*ndct[i] - common3[i]) -
                       m4p*common2[i]; 
            I_113[i] =  tmp[0]*bpctct[i];
            I_005[i] = -a2m8p*common2[i] - a2m4pn*common3[i] - tmp[1]*ndct[i];
            I_115[i] = -tmp[2]*bpctct[i] - tmp[3]*ndct[i];
        }
                     
        Fint_003 = f_103 - y[0]*f_003;
        Fint_113 = f_213 - y[0]*f_113;
        Fint_005 = f_105 - y[0]*f_005;
        Fint_115 = f_215 - y[0]*f_115;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            f4[i] = (I_003[i]*Fint_003 + I_113[i]*Fint_113 + I_005[i]*Fint_005 +
                     I_115[i]*Fint_115) * oneoverL;
        }

        Fint_003 = y[1]*f_003 - f_103;
        Fint_113 = y[1]*f_113 - f_213;
        Fint_005 = y[1]*f_005 - f_105;
        Fint_115 = y[1]*f_115 - f_215;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            f3[i] = (I_003[i]*Fint_003 + I_113[i]*Fint_113 + I_005[i]*Fint_005 +
                     I_115[i]*Fint_115) * oneoverL;
        }   
        
        fn3.x = f3[0];
        fn3.y = f3[1];
        fn3.z = f3[2];
        fn4.x = f4[0];
        fn4.y = f4[1];
        fn4.z = f4[2];
        
        x1modMag2 = 0.0e0;
        x2modMag2 = 0.0e0;
        
        #pragma unroll
        for (i=0;i<3;i++) {
            x1modMag2 += x1mod[i]*x1mod[i];
            x2modMag2 += x2mod[i]*x2mod[i];
        }

        diffMag2 = magdiff*magdiff;
        
        if (diffMag2 > (eps * (x1modMag2+x2modMag2))) {
			
			double3 fn3cor, fn4cor, fw, fq, bx1, bx2;
			double3 rx1, rx2, rx3, rx4, rx1mod, rx2mod;
			
			rx1.x = x1[0]; rx1.y = x1[1]; rx1.z = x1[2];
			rx2.x = x2[0]; rx2.y = x2[1]; rx2.z = x2[2];
			rx3.x = x3[0]; rx3.y = x3[1]; rx3.z = x3[2];
			rx4.x = x4[0]; rx4.y = x4[1]; rx4.z = x4[2];
			
			rx1mod.x = x1mod[0]; rx1mod.y = x1mod[1]; rx1mod.z = x1mod[2];
			rx2mod.x = x2mod[0]; rx2mod.y = x2mod[1]; rx2mod.z = x2mod[2];
			
			bx1.x = bp[0]; bx1.y = bp[1]; bx1.z = bp[2];
			bx2.x = b[0]; bx2.y = b[1]; bx2.z = b[2];
			
            SegSegForceIsotropicCorrGPU(rx1, rx1mod, rx3, rx4, bx1, bx2, fw, fq, fn3cor, fn4cor);
            
            fn3.x += fn3cor.x;
            fn3.y += fn3cor.y;
            fn3.z += fn3cor.z;
            fn4.x += fn4cor.x;
            fn4.y += fn4cor.y;
            fn4.z += fn4cor.z;
            
            SegSegForceIsotropicCorrGPU(rx2mod, rx2, rx3, rx4, bx1, bx2, fw, fq, fn3cor, fn4cor);
                       
            fn3.x += fn3cor.x;
            fn3.y += fn3cor.y;
            fn3.z += fn3cor.z;
            fn4.x += fn4cor.x;
            fn4.y += fn4cor.y;
            fn4.z += fn4cor.z;
        }
        
        return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    SegSegForceIsotropicCorrGPU
 *
 *-----------------------------------------------------------------------*/
__device__ void SegSegForceIsotropicCorrGPU(double3 r1, double3 r2, double3 r3, double3 r4, double3 b1, double3 b2,
                                            double3 &fn1, double3 &fn2, double3 &fn3, double3 &fn4)
{
        real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
        real8 f3[3], f4[3];
        real8 vec1[3], vec2[3], t[3], tp[3], tctp[3];
        real8 R[2][3], tempa[2], tempb[2], y[2], z[2];
        int i, j , alt1[3]={1,2,0}, alt2[3]={2,0,1};
        real8 d, c, c2, onemc2, onemc2inv, oneoverL, oneoverLp;
        real8 a2, m4p, m4pd, m8p, m8pd, m4pn, m4pnd, m4pnd2, m4pnd3;
        real8 a2m4pnd, a2m8pd, a2m4pn, a2m8p, a2_d2, a2_d2inv, denom;
        real8 temp1, temp2, temp3, temp4[8], tmp[10];
        real8 yv[4], zv[4], y2[4], z2[4], Ra[4], Rainv[4];
        real8 Ra_Rdot_tp[8], Ra_Rdot_t[8], log_Ra_Rdot_tp[4], log_Ra_Rdot_t[4];
        real8 Ra2_R_tpinv[4], Ra2_R_tinv[4], ylog_Ra_Rdot_tp[4], zlog_Ra_Rdot_t[4];
        real8 yRa2_R_tpinv[4], zRa2_R_tinv[4], y2Ra2_R_tpinv[4], z2Ra2_R_tinv[4];
        real8 adf_003[4], commonf223[4], commonf225[4], commonf025[4], commonf205[4];
        real8 commonf305[4], commonf035[4], ycommonf025[4], zcommonf205[4], zcommonf305[4];
        real8 tf_113[4];
        real8 f_003v[4], f_103v[4], f_013v[4], f_113v[4];
        real8 f_203v[4], f_023v[4], f_005v[4], f_105v[4];
        real8 f_003,  f_103,  f_013,  f_113,  f_203,  f_023,  f_005,  f_105;
        real8 f_015v[4], f_115v[4], f_205v[4], f_025v[4];
        real8 f_215v[4], f_125v[4], f_225v[4], f_305v[4];
        real8 f_015,  f_115,  f_205,  f_025,  f_215,  f_125,  f_225,  f_305;
        real8 f_035v[4], f_315v[4], f_135v[4];
        real8 f_035,  f_315,  f_135;
        real8 Fint_003, Fint_005, Fint_013, Fint_015, Fint_025, Fint_103;
        real8 Fint_105, Fint_115, Fint_125, Fint_205, Fint_215;
        real8 I_003[3], I_005[3], I_013[3], I_015[3], I_025[3], I_103[3];
        real8 I_105[3], I_115[3], I_125[3], I_205[3], I_215[3];
        real8 I00a[3], I01a[3], I10a[3], I00b[3], I01b[3], I10b[3];
        real8 bctctp[3], bct[3], bpctpct[3], bpctp[3], tcbpct[3];
        real8 bctdbp, bpctpdb, tcbpdb, tcbpdtp, tpcbdbp;
        real8 tctpct[3], tpct[3];
        real8 tctpcbpdb, tctpcbpdtp, tctpdb, tdb, tdbp;
        real8 tpcbctp[3], tpctctp[3];
        real8 tpcbdt, tpctcbdbp, tpctcbdt, tpctdbp, tpdb, tpdbp;
        real8 pivalue=3.141592653589793;           
        
        fn1.x = 0.0;
        fn1.y = 0.0;
        fn1.z = 0.0;

        fn2.x = 0.0;
        fn2.y = 0.0;
        fn2.z = 0.0;

        fn3.x = 0.0;
        fn3.y = 0.0;
        fn3.z = 0.0;

        fn4.x = 0.0;
        fn4.y = 0.0;
        fn4.z = 0.0;

        x1[0]=r1.x;
        x1[1]=r1.y;
        x1[2]=r1.z;
        x2[0]=r2.x;
        x2[1]=r2.y;
        x2[2]=r2.z;
        x3[0]=r3.x;
        x3[1]=r3.y;
        x3[2]=r3.z;
        x4[0]=r4.x;
        x4[1]=r4.y;
        x4[2]=r4.z;
        
        b[0]=b2.x;
        b[1]=b2.y;
        b[2]=b2.z;
        bp[0]=b1.x;
        bp[1]=b1.y;
        bp[2]=b1.z;
        
        #pragma unroll   
        for(i=0;i<3;i++) { 
            vec1[i]=x4[i]-x3[i];
            vec2[i]=x2[i]-x1[i];
        }

        temp1=0.0e0;
        temp2=0.0e0;    
        
        #pragma unroll 
        for(i=0;i<3;i++) { 
            temp1+=vec1[i]*vec1[i];
            temp2+=vec2[i]*vec2[i];
        }

        oneoverL =1/sqrt(temp1);
        oneoverLp=1/sqrt(temp2);
        
        #pragma unroll         
        for(i=0;i<3;i++) { 
            t[i]=vec1[i]*oneoverL;
            tp[i]=vec2[i]*oneoverLp;
        }
        
        c=0.0e0;
        
        #pragma unroll 
        for(i=0;i<3;i++) { 
            c+=t[i]*tp[i];
        }
                 
        c2=c*c;
        onemc2=1-c2;

        {	
			#pragma unroll 
            for(i=0;i<3;i++) {
                tctp[i]=t[alt1[i]]*tp[alt2[i]]-t[alt2[i]]*tp[alt1[i]];
            }

            onemc2inv = 1/onemc2;
            
            #pragma unroll 
            for(i=0;i<3;i++) { 
                R[0][i]=x3[i]-x1[i];
                R[1][i]=x4[i]-x2[i];
            }

            d=0.0e0;
            
            #pragma unroll
            for (j=0;j<2;j++) { 
                tempa[j]=0.0e0;
                tempb[j]=0.0e0;
            }
            
            #pragma unroll
            for(i=0;i<3;i++) { 
                d+=0.5e0*((x4[i]+x3[i])-(x2[i]+x1[i]))*tctp[i];
                for (j=0;j<2;j++) { 
                    tempa[j]+=R[j][i]*t[i];
                    tempb[j]+=R[j][i]*tp[i];
                }
            }

            d*=onemc2inv;
            
            #pragma unroll
            for (j=0;j<2;j++) { 
                y[j]=(tempa[j]-c*tempb[j])*onemc2inv;
                z[j]=(tempb[j]-c*tempa[j])*onemc2inv;
            }

/*          now we calculate the definite integrals of the force calculation  */

            #pragma unroll
            for (j=0;j<2;j++) {
                yv[2*j]=y[j];
                yv[2*j+1]=y[j];
                zv[j]=z[j];
                zv[j+2]=z[j];
            }
            
            a2_d2 = a*a+d*d*onemc2;
            
            #pragma unroll
            for (j=0;j<4;j++) {
                y2[j] = yv[j]*yv[j];
                z2[j] = zv[j]*zv[j];
                
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                temp4[j]=a2_d2 + y2[j] + z2[j] + 2.0e0*yv[j]*zv[j]*c;
            }

            temp1=onemc2*a2_d2;
            
            #pragma unroll
            for (j=0;j<4;j++) {
                Ra[j]=sqrt(temp4[j]);
            }

            temp2=sqrt(temp1);
            
            #pragma unroll
            for (j=0;j<4;j++) {
                Rainv[j]=1.0e0/Ra[j];
            }

            denom=1.0e0/temp2;
            a2_d2inv=1.0e0/a2_d2;
            
            #pragma unroll
            for (j=0;j<4;j++) {
                Ra_Rdot_tp[j] = Ra[j]+(zv[j]+yv[j]*c);       
                Ra_Rdot_t[j]  = Ra[j]+(yv[j]+zv[j]*c);
				Ra_Rdot_tp[j+4] = Ra[j]-(zv[j]+yv[j]*c);       
                Ra_Rdot_t[j+4]  = Ra[j]-(yv[j]+zv[j]*c); 
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                log_Ra_Rdot_tp[j] =0.5e0*(log(Ra_Rdot_tp[j])-log(Ra_Rdot_tp[j+4]));
                log_Ra_Rdot_t[j]  =0.5e0*(log(Ra_Rdot_t[j])-log(Ra_Rdot_t[j+4]));
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                Ra2_R_tpinv[j] = 0.5e0*(Rainv[j]/Ra_Rdot_tp[j]- Rainv[j]/Ra_Rdot_tp[j+4]);
                Ra2_R_tinv[j] =  0.5e0*(Rainv[j]/Ra_Rdot_t[j]- Rainv[j]/Ra_Rdot_t[j+4]);
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                ylog_Ra_Rdot_tp[j] = yv[j]*log_Ra_Rdot_tp[j];
                yRa2_R_tpinv[j]    = yv[j]*   Ra2_R_tpinv[j];
                zlog_Ra_Rdot_t[j]  = zv[j]*log_Ra_Rdot_t[j];
                zRa2_R_tinv[j]     = zv[j]*   Ra2_R_tinv[j];
                
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                y2Ra2_R_tpinv[j] = yv[j]* yRa2_R_tpinv[j];
                z2Ra2_R_tinv[j]  = zv[j]*  zRa2_R_tinv[j];
            }

            temp1=denom*(1+c);
            
            #pragma unroll
            for (j=0;j<4;j++) {
                temp4[j]=temp1*(Ra[j]+(yv[j]+zv[j]));
				temp4[j+4]=temp1*(Ra[j]-(yv[j]+zv[j]));
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                f_003v[j]=0.5e0*(atan(temp4[j])+atan(temp4[j+4]));
            }
            
            temp1=-2.0e0*denom;
            
            #pragma unroll
            for (j=0;j<4;j++) {
                f_003v[j]*=temp1;
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                adf_003[j]=f_003v[j]*a2_d2;
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                commonf223[j] = c*Ra[j] - adf_003[j];
                f_103v[j] = c*log_Ra_Rdot_t[j]  - log_Ra_Rdot_tp[j];
                f_013v[j] = c*log_Ra_Rdot_tp[j] - log_Ra_Rdot_t [j];
                f_113v[j] = c*adf_003[j] - Ra[j];
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                commonf223[j] *= onemc2inv;
                f_103v[j] *=      onemc2inv;
                f_013v[j] *=      onemc2inv;
                f_113v[j] *=      onemc2inv;
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                commonf225[j] = f_003v[j] - c*Rainv[j];
                commonf025[j] = c*yRa2_R_tpinv[j] - Rainv[j];
                commonf205[j] = c*zRa2_R_tinv[j]  - Rainv[j];
                commonf305[j] = log_Ra_Rdot_t[j]  -(yv[j]-c*zv[j])*Rainv[j] - c2*z2Ra2_R_tinv[j];
                commonf035[j] = log_Ra_Rdot_tp[j] -(zv[j]-c*yv[j])*Rainv[j] - c2*y2Ra2_R_tpinv[j]; 
                f_203v[j] =  zlog_Ra_Rdot_t[j]  + commonf223[j];
                f_023v[j] =  ylog_Ra_Rdot_tp[j] + commonf223[j];
                f_005v[j] = f_003v[j] - yRa2_R_tpinv[j] - zRa2_R_tinv[j];
                f_105v[j] = Ra2_R_tpinv[j] - c*Ra2_R_tinv[j];
                f_015v[j] = Ra2_R_tinv[j]  - c*Ra2_R_tpinv[j];
                f_115v[j] = Rainv[j] - c*(yRa2_R_tpinv[j] + zRa2_R_tinv[j] + f_003v[j]);
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                ycommonf025[j] = yv[j]*commonf025[j];
                zcommonf205[j] = zv[j]*commonf205[j];
                zcommonf305[j] = zv[j]*commonf305[j];
                tf_113[j]=2.0e0*f_113v[j];
                f_205v[j] = yRa2_R_tpinv[j] + c2*zRa2_R_tinv[j]  + commonf225[j];
                f_025v[j] = zRa2_R_tinv[j]  + c2*yRa2_R_tpinv[j] + commonf225[j];
                f_305v[j] = y2Ra2_R_tpinv[j] + c*commonf305[j] + 2.0e0*f_103v[j];
                f_035v[j] = z2Ra2_R_tinv[j]  + c*commonf035[j] + 2.0e0*f_013v[j];
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                f_215v[j] = f_013v[j] - ycommonf025[j] + c*(zcommonf205[j]-f_103v[j]); 
                f_125v[j] = f_103v[j] - zcommonf205[j] + c*(ycommonf025[j] - f_013v[j]); 
                f_225v[j] = f_203v[j] - zcommonf305[j] + c*(y2[j]*commonf025[j] - tf_113[j]);
                f_315v[j] = tf_113[j] - y2[j]*commonf025[j] + c*(zcommonf305[j] - f_203v[j]);
                f_135v[j] = tf_113[j] - z2[j]*commonf205[j] + c*(yv[j]*commonf035[j]-f_023v[j]);
            }
            
             
            f_003= (f_003v[0]+f_003v[3])-(f_003v[1]+f_003v[2]);
            f_013= (f_013v[0]+f_013v[3])-(f_013v[1]+f_013v[2]);
            f_103= (f_103v[0]+f_103v[3])-(f_103v[1]+f_103v[2]);
            f_113= (f_113v[0]+f_113v[3])-(f_113v[1]+f_113v[2]);
            f_023= (f_023v[0]+f_023v[3])-(f_023v[1]+f_023v[2]);
            f_203= (f_203v[0]+f_203v[3])-(f_203v[1]+f_203v[2]);
            f_005= (f_005v[0]+f_005v[3])-(f_005v[1]+f_005v[2]);
            f_015= (f_015v[0]+f_015v[3])-(f_015v[1]+f_015v[2]);
            f_105= (f_105v[0]+f_105v[3])-(f_105v[1]+f_105v[2]);
            f_115= (f_115v[0]+f_115v[3])-(f_115v[1]+f_115v[2]);
            f_025= (f_025v[0]+f_025v[3])-(f_025v[1]+f_025v[2]);
            f_205= (f_205v[0]+f_205v[3])-(f_205v[1]+f_205v[2]);
            f_215= (f_215v[0]+f_215v[3])-(f_215v[1]+f_215v[2]);
            f_125= (f_125v[0]+f_125v[3])-(f_125v[1]+f_125v[2]);
            f_035= (f_035v[0]+f_035v[3])-(f_035v[1]+f_035v[2]);
            f_305= (f_305v[0]+f_305v[3])-(f_305v[1]+f_305v[2]);
            f_225= (f_225v[0]+f_225v[3])-(f_225v[1]+f_225v[2]);
            f_135= (f_135v[0]+f_135v[3])-(f_135v[1]+f_135v[2]);
            f_315= (f_315v[0]+f_315v[3])-(f_315v[1]+f_315v[2]);
            
            
            f_005 *= a2_d2inv;
            f_105 *= onemc2inv;
            f_015 *= onemc2inv;
            f_115 *= onemc2inv;
            f_205 *= onemc2inv;
            f_025 *= onemc2inv;
            f_305 *= onemc2inv;
            f_035 *= onemc2inv;            
            f_215 *= onemc2inv; 
            f_125 *= onemc2inv; 
            f_225 *= onemc2inv;
            f_315 *= onemc2inv;
            f_135 *= onemc2inv;
            
      
/* now construct the vector coefficients for the definite integrals */

            a2 = a*a;
            m4p = 0.25 * MU / pivalue;
            m4pd =  m4p * d;
            m8p = 0.5 * m4p;
            m8pd = m8p * d;
            m4pn = m4p / ( 1 - NU );
            m4pnd = m4pn * d;
            m4pnd2 = m4pnd * d;
            m4pnd3 = m4pnd2 * d;
            a2m4pnd = a2 * m4pnd;
            a2m8pd = a2 * m8pd;
            a2m4pn = a2 * m4pn;
            a2m8p = a2 * m8p;

            #pragma unroll
            for (i=0;i<3;i++) {
                tpct[i]=-tctp[i];
                bct[i]=b[alt1[i]]*t[alt2[i]]-b[alt2[i]]*t[alt1[i]];
                bpctp[i]=bp[alt1[i]]*tp[alt2[i]]-bp[alt2[i]]*tp[alt1[i]];
                
            }

            tdb=0.0e0;
            tdbp=0.0e0;
            tpdb=0.0e0;
            tpdbp=0.0e0;
            tctpdb=0.0e0;
            tpctdbp=0.0e0;
            bpctpdb=0.0e0;
            bctdbp=0.0e0;
            
            #pragma unroll
            for (i=0;i<3;i++) {
                tdb    +=t[i]*b[i];
                tdbp   +=t[i]*bp[i];
                tpdb   +=tp[i]*b[i];
                tpdbp  +=tp[i]*bp[i];
                tctpdb +=tctp[i]*b[i];
                tpctdbp+=tpct[i]*bp[i];
                bpctpdb+=bpctp[i]*b[i];
                bctdbp +=bct[i]*bp[i];
            }
            
            #pragma unroll
            for (i=0;i<3;i++) {
                tctpct[i]    =        tp[i] -     c*t[i];
                tpctctp[i]   =         t[i] -    c*tp[i];
                tcbpct[i]    =        bp[i] -  tdbp*t[i];
                tpcbctp[i]   =         b[i] - tpdb*tp[i];
                bpctpct[i]   =   tdbp*tp[i] -    c*bp[i];
                bctctp[i]    =    tpdb*t[i] -     c*b[i];
            }
                
            
            tctpcbpdtp = tdbp - tpdbp*c;
            tpctcbdt = tpdb - tdb*c;
            tctpcbpdb =  tdbp*tpdb - tpdbp*tdb;
            tpctcbdbp = tctpcbpdb;
            tcbpdtp = tpctdbp; 
            tpcbdt = tctpdb;
            tcbpdb = bctdbp;
            tpcbdbp = bpctpdb;

/*
 *          Only calculate the forces for segment p3->p4 if at least one
 *          of the segment's nodes is local to the current domain.
 */
            //if (seg34Local) {

                temp1 = tdbp*tpdb + tctpcbpdb;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I00a[i] = temp1 * tpct[i];
                    I00b[i] = tctpcbpdtp * bct[i];
                }

                temp1 = (m4pnd * tctpdb);
                temp2 = (m4pnd * bpctpdb);
                temp3 = (m4pnd3 * tctpcbpdtp*tctpdb);
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_003[i] = m4pd*I00a[i] - m4pnd*I00b[i] + temp1*bpctpct[i] +
                            temp2*tctpct[i]; 
                    I_005[i] = a2m8pd*I00a[i] - a2m4pnd*I00b[i] - temp3*tctpct[i];
                    I10a[i] = tcbpct[i]*tpdb - tctp[i]*tcbpdb;
                    I10b[i] = bct[i] * tcbpdtp;
                    
                }

                temp1 = (m4pn * tdb);
                temp2 = m4pnd2 * (tcbpdtp*tctpdb + tctpcbpdtp*tdb);
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_103[i] = temp1*bpctpct[i] + m4p*I10a[i] - m4pn*I10b[i];
                    I_105[i] = a2m8p*I10a[i] - a2m4pn*I10b[i] - temp2*tctpct[i];
                    I01a[i] = tctp[i]*bpctpdb - bpctpct[i]*tpdb;
                }

                tmp[0] = (m4pn * tpdb); 
                tmp[1] = (m4pn * bpctpdb);
                tmp[2] = (m4pnd2 * tctpcbpdtp * tpdb);
                tmp[3] = (m4pnd2 * tctpcbpdtp * tctpdb);
                tmp[4] = (m4pnd * tcbpdtp * tdb);
                tmp[5] = (m4pnd * tctpcbpdtp * tpdb) ;
                tmp[6] = (m4pnd * (tctpcbpdtp*tdb + tcbpdtp*tctpdb));
                tmp[7] = (m4pnd * tcbpdtp * tpdb);
                tmp[8] = (m4pn * tcbpdtp * tdb);
                tmp[9] = (m4pn * tcbpdtp * tpdb);
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_013[i] = m4p*I01a[i] + tmp[0]*bpctpct[i] - tmp[1]*tctp[i];
                    I_015[i] = a2m8p*I01a[i] - tmp[2]*tctpct[i] + tmp[3]*tctp[i];
                    I_205[i] = -tmp[4] * tctpct[i];
                    I_025[i] = tmp[5] * tctp[i]; 
                    I_115[i] =  tmp[6]*tctp[i] - tmp[7]*tctpct[i];
                    I_215[i] = tmp[8] * tctp[i];
                    I_125[i] = tmp[9] * tctp[i];
                }
  
                Fint_003 = f_103 - y[0]*f_003;
                Fint_103 = f_203 - y[0]*f_103;
                Fint_013 = f_113 - y[0]*f_013;
                Fint_005 = f_105 - y[0]*f_005;
                Fint_105 = f_205 - y[0]*f_105;
                Fint_015 = f_115 - y[0]*f_015;
                Fint_115 = f_215 - y[0]*f_115;
                Fint_205 = f_305 - y[0]*f_205;
                Fint_025 = f_125 - y[0]*f_025;
                Fint_215 = f_315 - y[0]*f_215;
                Fint_125 = f_225 - y[0]*f_125;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    f4[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverL;
                }

                Fint_003 = y[1]*f_003 - f_103;
                Fint_103 = y[1]*f_103 - f_203;
                Fint_013 = y[1]*f_013 - f_113;
                Fint_005 = y[1]*f_005 - f_105;
                Fint_105 = y[1]*f_105 - f_205;
                Fint_015 = y[1]*f_015 - f_115;
                Fint_115 = y[1]*f_115 - f_215;
                Fint_205 = y[1]*f_205 - f_305;
                Fint_025 = y[1]*f_025 - f_125;
                Fint_215 = y[1]*f_215 - f_315;
                Fint_125 = y[1]*f_125 - f_225;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    f3[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverL;
                }

                fn3.x=f3[0];
                fn3.y=f3[1];
                fn3.z=f3[2];
                fn4.x=f4[0];
                fn4.y=f4[1];
                fn4.z=f4[2];

            //} /* if segment p3->p4 is "local" */
       }
       
       return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    SegSegForceIsotropicGPU
 *
 *-----------------------------------------------------------------------*/
__device__ void SegSegForceIsotropicGPU(double3 r1, double3 r2, double3 r3, double3 r4, double3 b1, double3 b2,
                                        double3 &fn1, double3 &fn2, double3 &fn3, double3 &fn4)
{
        real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
        real8 f1[3], f2[3], f3[3], f4[3];
        real8 vec1[3], vec2[3], t[3], tp[3], tctp[3];
        real8 R[2][3], tempa[2], tempb[2], y[2], z[2];
        int i, j , alt1[3]={1,2,0}, alt2[3]={2,0,1};
        real8 eps, d, c, c2, onemc2, onemc2inv, oneoverL, oneoverLp;
        real8 a2, m4p, m4pd, m8p, m8pd, m4pn, m4pnd, m4pnd2, m4pnd3;
        real8 a2m4pnd, a2m8pd, a2m4pn, a2m8p, a2_d2, a2_d2inv, denom;
        real8 temp1, temp2, temp3, temp4[8], tmp[10];
        real8 yv[4], zv[4], y2[4], z2[4], Ra[4], Rainv[4];
        real8 Ra_Rdot_tp[8], Ra_Rdot_t[8], log_Ra_Rdot_tp[4], log_Ra_Rdot_t[4];
        real8 Ra2_R_tpinv[4], Ra2_R_tinv[4], ylog_Ra_Rdot_tp[4], zlog_Ra_Rdot_t[4];
        real8 yRa2_R_tpinv[4], zRa2_R_tinv[4], y2Ra2_R_tpinv[4], z2Ra2_R_tinv[4];
        real8 adf_003[4], commonf223[4], commonf225[4], commonf025[4], commonf205[4];
        real8 commonf305[4], commonf035[4], ycommonf025[4], zcommonf205[4], zcommonf305[4];
        real8 tf_113[4];
        real8 f_003v[4], f_103v[4], f_013v[4], f_113v[4];
        real8 f_203v[4], f_023v[4], f_005v[4], f_105v[4];
        real8 f_003,  f_103,  f_013,  f_113,  f_203,  f_023,  f_005,  f_105;
        real8 f_015v[4], f_115v[4], f_205v[4], f_025v[4];
        real8 f_215v[4], f_125v[4], f_225v[4], f_305v[4];
        real8 f_015,  f_115,  f_205,  f_025,  f_215,  f_125,  f_225,  f_305;
        real8 f_035v[4], f_315v[4], f_135v[4];
        real8 f_035,  f_315,  f_135;
        real8 Fint_003, Fint_005, Fint_013, Fint_015, Fint_025, Fint_103;
        real8 Fint_105, Fint_115, Fint_125, Fint_205, Fint_215;
        real8 I_003[3], I_005[3], I_013[3], I_015[3], I_025[3], I_103[3];
        real8 I_105[3], I_115[3], I_125[3], I_205[3], I_215[3];
        real8 I00a[3], I01a[3], I10a[3], I00b[3], I01b[3], I10b[3];
        real8 bctctp[3], bct[3], bpctpct[3], bpctp[3], tcbpct[3];
        real8 bctdbp, bpctpdb, tcbpdb, tcbpdtp, tpcbdbp;
        real8 tctpct[3], tpct[3];
        real8 tctpcbpdb, tctpcbpdtp, tctpdb, tdb, tdbp;
        real8 tpcbctp[3], tpctctp[3];
        real8 tpcbdt, tpctcbdbp, tpctcbdt, tpctdbp, tpdb, tpdbp;
        real8 pivalue=3.141592653589793;

        eps = 1e-4;            
        
        fn1.x = 0.0;
        fn1.y = 0.0;
        fn1.z = 0.0;

        fn2.x = 0.0;
        fn2.y = 0.0;
        fn2.z = 0.0;

        fn3.x = 0.0;
        fn3.y = 0.0;
        fn3.z = 0.0;

        fn4.x = 0.0;
        fn4.y = 0.0;
        fn4.z = 0.0;

        x1[0]=r1.x;
        x1[1]=r1.y;
        x1[2]=r1.z;
        x2[0]=r2.x;
        x2[1]=r2.y;
        x2[2]=r2.z;
        x3[0]=r3.x;
        x3[1]=r3.y;
        x3[2]=r3.z;
        x4[0]=r4.x;
        x4[1]=r4.y;
        x4[2]=r4.z;
        
        b[0]=b2.x;
        b[1]=b2.y;
        b[2]=b2.z;
        bp[0]=b1.x;
        bp[1]=b1.y;
        bp[2]=b1.z;
        
        #pragma unroll   
        for(i=0;i<3;i++) { 
            vec1[i]=x4[i]-x3[i];
            vec2[i]=x2[i]-x1[i];
        }

        temp1=0.0e0;
        temp2=0.0e0;    
        
        #pragma unroll 
        for(i=0;i<3;i++) { 
            temp1+=vec1[i]*vec1[i];
            temp2+=vec2[i]*vec2[i];
        }

        oneoverL =1/sqrt(temp1);
        oneoverLp=1/sqrt(temp2);
        
        #pragma unroll         
        for(i=0;i<3;i++) { 
            t[i]=vec1[i]*oneoverL;
            tp[i]=vec2[i]*oneoverLp;
        }
        
        c=0.0e0;
        
        #pragma unroll 
        for(i=0;i<3;i++) { 
            c+=t[i]*tp[i];
        }
                 
        c2=c*c;
        onemc2=1-c2;
        
        if (onemc2 > eps) {
			
			#pragma unroll 
            for(i=0;i<3;i++) {
                tctp[i]=t[alt1[i]]*tp[alt2[i]]-t[alt2[i]]*tp[alt1[i]];
            }

            onemc2inv = 1/onemc2;
            
            #pragma unroll 
            for(i=0;i<3;i++) { 
                R[0][i]=x3[i]-x1[i];
                R[1][i]=x4[i]-x2[i];
            }

            d=0.0e0;
            
            #pragma unroll
            for (j=0;j<2;j++) { 
                tempa[j]=0.0e0;
                tempb[j]=0.0e0;
            }
            
            #pragma unroll
            for(i=0;i<3;i++) { 
                d+=0.5e0*((x4[i]+x3[i])-(x2[i]+x1[i]))*tctp[i];
                for (j=0;j<2;j++) { 
                    tempa[j]+=R[j][i]*t[i];
                    tempb[j]+=R[j][i]*tp[i];
                }
            }

            d*=onemc2inv;
            
            #pragma unroll
            for (j=0;j<2;j++) { 
                y[j]=(tempa[j]-c*tempb[j])*onemc2inv;
                z[j]=(tempb[j]-c*tempa[j])*onemc2inv;
            }

/*          now we calculate the definite integrals of the force calculation  */

            #pragma unroll
            for (j=0;j<2;j++) {
                yv[2*j]=y[j];
                yv[2*j+1]=y[j];
                zv[j]=z[j];
                zv[j+2]=z[j];
            }
            
            a2_d2 = a*a+d*d*onemc2;
            
            #pragma unroll
            for (j=0;j<4;j++) {
                y2[j] = yv[j]*yv[j];
                z2[j] = zv[j]*zv[j];
                
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                temp4[j]=a2_d2 + y2[j] + z2[j] + 2.0e0*yv[j]*zv[j]*c;
            }

            temp1=onemc2*a2_d2;
            
            #pragma unroll
            for (j=0;j<4;j++) {
                Ra[j]=sqrt(temp4[j]);
            }

            temp2=sqrt(temp1);
            
            #pragma unroll
            for (j=0;j<4;j++) {
                Rainv[j]=1.0e0/Ra[j];
            }

            denom=1.0e0/temp2;
            a2_d2inv=1.0e0/a2_d2;
            
            #pragma unroll
            for (j=0;j<4;j++) {
                Ra_Rdot_tp[j] = Ra[j]+(zv[j]+yv[j]*c);       
                Ra_Rdot_t[j]  = Ra[j]+(yv[j]+zv[j]*c);
				Ra_Rdot_tp[j+4] = Ra[j]-(zv[j]+yv[j]*c);       
                Ra_Rdot_t[j+4]  = Ra[j]-(yv[j]+zv[j]*c); 
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                log_Ra_Rdot_tp[j] =0.5e0*(log(Ra_Rdot_tp[j])-log(Ra_Rdot_tp[j+4]));
                log_Ra_Rdot_t[j]  =0.5e0*(log(Ra_Rdot_t[j])-log(Ra_Rdot_t[j+4]));
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                Ra2_R_tpinv[j] = 0.5e0*(Rainv[j]/Ra_Rdot_tp[j]- Rainv[j]/Ra_Rdot_tp[j+4]);
                Ra2_R_tinv[j] =  0.5e0*(Rainv[j]/Ra_Rdot_t[j]- Rainv[j]/Ra_Rdot_t[j+4]);
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                ylog_Ra_Rdot_tp[j] = yv[j]*log_Ra_Rdot_tp[j];
                yRa2_R_tpinv[j]    = yv[j]*   Ra2_R_tpinv[j];
                zlog_Ra_Rdot_t[j]  = zv[j]*log_Ra_Rdot_t[j];
                zRa2_R_tinv[j]     = zv[j]*   Ra2_R_tinv[j];
                
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                y2Ra2_R_tpinv[j] = yv[j]* yRa2_R_tpinv[j];
                z2Ra2_R_tinv[j]  = zv[j]*  zRa2_R_tinv[j];
            }

            temp1=denom*(1+c);
            
            #pragma unroll
            for (j=0;j<4;j++) {
                temp4[j]=temp1*(Ra[j]+(yv[j]+zv[j]));
				temp4[j+4]=temp1*(Ra[j]-(yv[j]+zv[j]));
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                f_003v[j]=0.5e0*(atan(temp4[j])+atan(temp4[j+4]));
            }
            
            temp1=-2.0e0*denom;
            
            #pragma unroll
            for (j=0;j<4;j++) {
                f_003v[j]*=temp1;
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                adf_003[j]=f_003v[j]*a2_d2;
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                commonf223[j] = c*Ra[j] - adf_003[j];
                f_103v[j] = c*log_Ra_Rdot_t[j]  - log_Ra_Rdot_tp[j];
                f_013v[j] = c*log_Ra_Rdot_tp[j] - log_Ra_Rdot_t [j];
                f_113v[j] = c*adf_003[j] - Ra[j];
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                commonf223[j] *= onemc2inv;
                f_103v[j] *=      onemc2inv;
                f_013v[j] *=      onemc2inv;
                f_113v[j] *=      onemc2inv;
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                commonf225[j] = f_003v[j] - c*Rainv[j];
                commonf025[j] = c*yRa2_R_tpinv[j] - Rainv[j];
                commonf205[j] = c*zRa2_R_tinv[j]  - Rainv[j];
                commonf305[j] = log_Ra_Rdot_t[j]  -(yv[j]-c*zv[j])*Rainv[j] - c2*z2Ra2_R_tinv[j];
                commonf035[j] = log_Ra_Rdot_tp[j] -(zv[j]-c*yv[j])*Rainv[j] - c2*y2Ra2_R_tpinv[j]; 
                f_203v[j] =  zlog_Ra_Rdot_t[j]  + commonf223[j];
                f_023v[j] =  ylog_Ra_Rdot_tp[j] + commonf223[j];
                f_005v[j] = f_003v[j] - yRa2_R_tpinv[j] - zRa2_R_tinv[j];
                f_105v[j] = Ra2_R_tpinv[j] - c*Ra2_R_tinv[j];
                f_015v[j] = Ra2_R_tinv[j]  - c*Ra2_R_tpinv[j];
                f_115v[j] = Rainv[j] - c*(yRa2_R_tpinv[j] + zRa2_R_tinv[j] + f_003v[j]);
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                ycommonf025[j] = yv[j]*commonf025[j];
                zcommonf205[j] = zv[j]*commonf205[j];
                zcommonf305[j] = zv[j]*commonf305[j];
                tf_113[j]=2.0e0*f_113v[j];
                f_205v[j] = yRa2_R_tpinv[j] + c2*zRa2_R_tinv[j]  + commonf225[j];
                f_025v[j] = zRa2_R_tinv[j]  + c2*yRa2_R_tpinv[j] + commonf225[j];
                f_305v[j] = y2Ra2_R_tpinv[j] + c*commonf305[j] + 2.0e0*f_103v[j];
                f_035v[j] = z2Ra2_R_tinv[j]  + c*commonf035[j] + 2.0e0*f_013v[j];
            }
            
            #pragma unroll
            for (j=0;j<4;j++) {
                f_215v[j] = f_013v[j] - ycommonf025[j] + c*(zcommonf205[j]-f_103v[j]); 
                f_125v[j] = f_103v[j] - zcommonf205[j] + c*(ycommonf025[j] - f_013v[j]); 
                f_225v[j] = f_203v[j] - zcommonf305[j] + c*(y2[j]*commonf025[j] - tf_113[j]);
                f_315v[j] = tf_113[j] - y2[j]*commonf025[j] + c*(zcommonf305[j] - f_203v[j]);
                f_135v[j] = tf_113[j] - z2[j]*commonf205[j] + c*(yv[j]*commonf035[j]-f_023v[j]);
            }
            
             
            f_003= (f_003v[0]+f_003v[3])-(f_003v[1]+f_003v[2]);
            f_013= (f_013v[0]+f_013v[3])-(f_013v[1]+f_013v[2]);
            f_103= (f_103v[0]+f_103v[3])-(f_103v[1]+f_103v[2]);
            f_113= (f_113v[0]+f_113v[3])-(f_113v[1]+f_113v[2]);
            f_023= (f_023v[0]+f_023v[3])-(f_023v[1]+f_023v[2]);
            f_203= (f_203v[0]+f_203v[3])-(f_203v[1]+f_203v[2]);
            f_005= (f_005v[0]+f_005v[3])-(f_005v[1]+f_005v[2]);
            f_015= (f_015v[0]+f_015v[3])-(f_015v[1]+f_015v[2]);
            f_105= (f_105v[0]+f_105v[3])-(f_105v[1]+f_105v[2]);
            f_115= (f_115v[0]+f_115v[3])-(f_115v[1]+f_115v[2]);
            f_025= (f_025v[0]+f_025v[3])-(f_025v[1]+f_025v[2]);
            f_205= (f_205v[0]+f_205v[3])-(f_205v[1]+f_205v[2]);
            f_215= (f_215v[0]+f_215v[3])-(f_215v[1]+f_215v[2]);
            f_125= (f_125v[0]+f_125v[3])-(f_125v[1]+f_125v[2]);
            f_035= (f_035v[0]+f_035v[3])-(f_035v[1]+f_035v[2]);
            f_305= (f_305v[0]+f_305v[3])-(f_305v[1]+f_305v[2]);
            f_225= (f_225v[0]+f_225v[3])-(f_225v[1]+f_225v[2]);
            f_135= (f_135v[0]+f_135v[3])-(f_135v[1]+f_135v[2]);
            f_315= (f_315v[0]+f_315v[3])-(f_315v[1]+f_315v[2]);
            
            
            f_005 *= a2_d2inv;
            f_105 *= onemc2inv;
            f_015 *= onemc2inv;
            f_115 *= onemc2inv;
            f_205 *= onemc2inv;
            f_025 *= onemc2inv;
            f_305 *= onemc2inv;
            f_035 *= onemc2inv;            
            f_215 *= onemc2inv; 
            f_125 *= onemc2inv; 
            f_225 *= onemc2inv;
            f_315 *= onemc2inv;
            f_135 *= onemc2inv;
            
      
/* now construct the vector coefficients for the definite integrals */

            a2 = a*a;
            m4p = 0.25 * MU / pivalue;
            m4pd =  m4p * d;
            m8p = 0.5 * m4p;
            m8pd = m8p * d;
            m4pn = m4p / ( 1 - NU );
            m4pnd = m4pn * d;
            m4pnd2 = m4pnd * d;
            m4pnd3 = m4pnd2 * d;
            a2m4pnd = a2 * m4pnd;
            a2m8pd = a2 * m8pd;
            a2m4pn = a2 * m4pn;
            a2m8p = a2 * m8p;

            #pragma unroll
            for (i=0;i<3;i++) {
                tpct[i]=-tctp[i];
                bct[i]=b[alt1[i]]*t[alt2[i]]-b[alt2[i]]*t[alt1[i]];
                bpctp[i]=bp[alt1[i]]*tp[alt2[i]]-bp[alt2[i]]*tp[alt1[i]];
                
            }

            tdb=0.0e0;
            tdbp=0.0e0;
            tpdb=0.0e0;
            tpdbp=0.0e0;
            tctpdb=0.0e0;
            tpctdbp=0.0e0;
            bpctpdb=0.0e0;
            bctdbp=0.0e0;
            
            #pragma unroll
            for (i=0;i<3;i++) {
                tdb    +=t[i]*b[i];
                tdbp   +=t[i]*bp[i];
                tpdb   +=tp[i]*b[i];
                tpdbp  +=tp[i]*bp[i];
                tctpdb +=tctp[i]*b[i];
                tpctdbp+=tpct[i]*bp[i];
                bpctpdb+=bpctp[i]*b[i];
                bctdbp +=bct[i]*bp[i];
            }
            
            #pragma unroll
            for (i=0;i<3;i++) {
                tctpct[i]    =        tp[i] -     c*t[i];
                tpctctp[i]   =         t[i] -    c*tp[i];
                tcbpct[i]    =        bp[i] -  tdbp*t[i];
                tpcbctp[i]   =         b[i] - tpdb*tp[i];
                bpctpct[i]   =   tdbp*tp[i] -    c*bp[i];
                bctctp[i]    =    tpdb*t[i] -     c*b[i];
            }
                
            
            tctpcbpdtp = tdbp - tpdbp*c;
            tpctcbdt = tpdb - tdb*c;
            tctpcbpdb =  tdbp*tpdb - tpdbp*tdb;
            tpctcbdbp = tctpcbpdb;
            tcbpdtp = tpctdbp; 
            tpcbdt = tctpdb;
            tcbpdb = bctdbp;
            tpcbdbp = bpctpdb;

/*
 *          Only calculate the forces for segment p3->p4 if at least one
 *          of the segment's nodes is local to the current domain.
 */
            //if (seg34Local) {

                temp1 = tdbp*tpdb + tctpcbpdb;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I00a[i] = temp1 * tpct[i];
                    I00b[i] = tctpcbpdtp * bct[i];
                }

                temp1 = (m4pnd * tctpdb);
                temp2 = (m4pnd * bpctpdb);
                temp3 = (m4pnd3 * tctpcbpdtp*tctpdb);
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_003[i] = m4pd*I00a[i] - m4pnd*I00b[i] + temp1*bpctpct[i] +
                            temp2*tctpct[i]; 
                    I_005[i] = a2m8pd*I00a[i] - a2m4pnd*I00b[i] - temp3*tctpct[i];
                    I10a[i] = tcbpct[i]*tpdb - tctp[i]*tcbpdb;
                    I10b[i] = bct[i] * tcbpdtp;
                    
                }

                temp1 = (m4pn * tdb);
                temp2 = m4pnd2 * (tcbpdtp*tctpdb + tctpcbpdtp*tdb);
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_103[i] = temp1*bpctpct[i] + m4p*I10a[i] - m4pn*I10b[i];
                    I_105[i] = a2m8p*I10a[i] - a2m4pn*I10b[i] - temp2*tctpct[i];
                    I01a[i] = tctp[i]*bpctpdb - bpctpct[i]*tpdb;
                }

                tmp[0] = (m4pn * tpdb); 
                tmp[1] = (m4pn * bpctpdb);
                tmp[2] = (m4pnd2 * tctpcbpdtp * tpdb);
                tmp[3] = (m4pnd2 * tctpcbpdtp * tctpdb);
                tmp[4] = (m4pnd * tcbpdtp * tdb);
                tmp[5] = (m4pnd * tctpcbpdtp * tpdb) ;
                tmp[6] = (m4pnd * (tctpcbpdtp*tdb + tcbpdtp*tctpdb));
                tmp[7] = (m4pnd * tcbpdtp * tpdb);
                tmp[8] = (m4pn * tcbpdtp * tdb);
                tmp[9] = (m4pn * tcbpdtp * tpdb);
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_013[i] = m4p*I01a[i] + tmp[0]*bpctpct[i] - tmp[1]*tctp[i];
                    I_015[i] = a2m8p*I01a[i] - tmp[2]*tctpct[i] + tmp[3]*tctp[i];
                    I_205[i] = -tmp[4] * tctpct[i];
                    I_025[i] = tmp[5] * tctp[i]; 
                    I_115[i] =  tmp[6]*tctp[i] - tmp[7]*tctpct[i];
                    I_215[i] = tmp[8] * tctp[i];
                    I_125[i] = tmp[9] * tctp[i];
                }
  
                Fint_003 = f_103 - y[0]*f_003;
                Fint_103 = f_203 - y[0]*f_103;
                Fint_013 = f_113 - y[0]*f_013;
                Fint_005 = f_105 - y[0]*f_005;
                Fint_105 = f_205 - y[0]*f_105;
                Fint_015 = f_115 - y[0]*f_015;
                Fint_115 = f_215 - y[0]*f_115;
                Fint_205 = f_305 - y[0]*f_205;
                Fint_025 = f_125 - y[0]*f_025;
                Fint_215 = f_315 - y[0]*f_215;
                Fint_125 = f_225 - y[0]*f_125;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    f4[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverL;
                }

                Fint_003 = y[1]*f_003 - f_103;
                Fint_103 = y[1]*f_103 - f_203;
                Fint_013 = y[1]*f_013 - f_113;
                Fint_005 = y[1]*f_005 - f_105;
                Fint_105 = y[1]*f_105 - f_205;
                Fint_015 = y[1]*f_015 - f_115;
                Fint_115 = y[1]*f_115 - f_215;
                Fint_205 = y[1]*f_205 - f_305;
                Fint_025 = y[1]*f_025 - f_125;
                Fint_215 = y[1]*f_215 - f_315;
                Fint_125 = y[1]*f_125 - f_225;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    f3[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverL;
                }

                fn3.x=f3[0];
                fn3.y=f3[1];
                fn3.z=f3[2];
                fn4.x=f4[0];
                fn4.y=f4[1];
                fn4.z=f4[2];

            //} /* if segment p3->p4 is "local" */

/*
 *          Only calculate the forces for segment p1->p2 if at least one
 *          of the segment's nodes is local to the current domain.
 */
            //if (seg12Local) {

                temp1 = tpdb*tdbp + tpctcbdbp;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I00a[i] = temp1 * tctp[i];
                    I00b[i] = bpctp[i] * tpctcbdt;
                }
                
                temp1 = m4pnd * tpctdbp;
                temp2 = m4pnd * bctdbp;
                temp3 = m4pnd3 * tpctcbdt * tpctdbp;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_003[i] = m4pd*I00a[i] - m4pnd*I00b[i] + temp1*bctctp[i] +
                               temp2*tpctctp[i];
                    I_005[i] = a2m8pd*I00a[i] - a2m4pnd*I00b[i] - temp3*tpctctp[i]; 
                    I01a[i] = tpct[i]*tpcbdbp - tpcbctp[i]*tdbp;
                    I01b[i] = -bpctp[i] * tpcbdt;
                }

                temp1 = m4pn * tpdbp;
                temp2 = m4pnd2 * (tpcbdt*tpctdbp + tpctcbdt*tpdbp);
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_013[i] = -temp1 * bctctp[i] + m4p*I01a[i] - m4pn*I01b[i];
                    I_015[i] = a2m8p*I01a[i] - a2m4pn*I01b[i] + temp2*tpctctp[i];
                    I10a[i] = bctctp[i]*tdbp - tpct[i]*bctdbp;
                }

                tmp[0] = m4pn * tdbp; 
                tmp[1] = m4pn * bctdbp;
                tmp[2] = m4pnd2 * tpctcbdt * tdbp;
                tmp[3] = m4pnd2 * tpctcbdt * tpctdbp;
                tmp[4] = (m4pnd * tpcbdt * tpdbp);
                tmp[5] = (m4pnd * tpctcbdt * tdbp);
                tmp[6] = m4pnd * (tpctcbdt*tpdbp + tpcbdt*tpctdbp);
                tmp[7] = m4pnd * tpcbdt * tdbp;
                tmp[8] = (m4pn * tpcbdt * tpdbp);
                tmp[9] = (m4pn * tpcbdt * tdbp);
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    I_103[i] = m4p*I10a[i] - tmp[0]*bctctp[i] + tmp[1]*tpct[i];
                    I_105[i] = a2m8p*I10a[i] + tmp[2]*tpctctp[i] - tmp[3]*tpct[i];
                    I_025[i] = -tmp[4] * tpctctp[i];
                    I_205[i] = tmp[5] * tpct[i];
                    I_115[i] = tmp[6]*tpct[i] - tmp[7]*tpctctp[i];
                    I_125[i] = -tmp[8] * tpct[i];
                    I_215[i] = -tmp[9] * tpct[i];
                }

                Fint_003 = f_013 - z[1]*f_003;
                Fint_103 = f_113 - z[1]*f_103;
                Fint_013 = f_023 - z[1]*f_013;
                Fint_005 = f_015 - z[1]*f_005;
                Fint_105 = f_115 - z[1]*f_105;
                Fint_015 = f_025 - z[1]*f_015;
                Fint_115 = f_125 - z[1]*f_115;
                Fint_205 = f_215 - z[1]*f_205;
                Fint_025 = f_035 - z[1]*f_025;
                Fint_215 = f_225 - z[1]*f_215;
                Fint_125 = f_135 - z[1]*f_125;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    f1[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverLp;
                }
   
                Fint_003 = z[0]*f_003 - f_013;
                Fint_103 = z[0]*f_103 - f_113;
                Fint_013 = z[0]*f_013 - f_023;
                Fint_005 = z[0]*f_005 - f_015;
                Fint_105 = z[0]*f_105 - f_115;
                Fint_015 = z[0]*f_015 - f_025;
                Fint_115 = z[0]*f_115 - f_125;
                Fint_205 = z[0]*f_205 - f_215;
                Fint_025 = z[0]*f_025 - f_035;
                Fint_215 = z[0]*f_215 - f_225;
                Fint_125 = z[0]*f_125 - f_135;
                
                #pragma unroll
                for (i=0;i<3;i++) {
                    f2[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 + I_013[i]*Fint_013 +
                           I_005[i]*Fint_005 + I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                           I_115[i]*Fint_115 + I_205[i]*Fint_205 + I_025[i]*Fint_025 +
                           I_215[i]*Fint_215 + I_125[i]*Fint_125) * oneoverLp;
                }
                
                fn1.x=f1[0];
                fn1.y=f1[1];
                fn1.z=f1[2];
                fn2.x=f2[0];
                fn2.y=f2[1];
                fn2.z=f2[2];
                
   
            //} /* if segment p1->p2 is "local" */

        } else {
/*
 *          The two lines are parallel, so we have to use a special
 *          lower dimensional function
 */
			
			SpecialSegSegForceHalfGPU(r1, r2, r3, r4, b1, b2, fn3, fn4);

            SpecialSegSegForceHalfGPU(r3, r4, r1, r2, b2, b1, fn1, fn2);
       }

       return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    DotProductGPU
 *
 *-----------------------------------------------------------------------*/
__device__ double DotProductGPU(double3 v1, double3 v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

/*------------------------------------------------------------------------
 *
 *      Function:    ZImageGPU
 *
 *-----------------------------------------------------------------------*/
__device__ void ZImageGPU(double3 &r)
{
/*
 *      If periodic boundaries are not in use, the provided position
 *      of (x,y,z) will not be adjusted since there are no other
 *      images available (in that case Lx,Ly,Lz = 0).
 */
		r.x -= rint(r.x * invLx) * Lx;
		r.y -= rint(r.y * invLy) * Ly;
		r.z -= rint(r.z * invLz) * Lz;
}

/*------------------------------------------------------------------------
 *
 *      Function:    FoldBoxGPU
 *
 *-----------------------------------------------------------------------*/
__device__ void FoldBoxGPU(double3 &r)
{
		r.x -= rint((r.x-boxc.x)*invLx) * Lx;
		r.y -= rint((r.y-boxc.y)*invLy) * Ly;
		r.z -= rint((r.z-boxc.z)*invLz) * Lz;
}

/*------------------------------------------------------------------------
 *
 *      Function:    PreserveNodesGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void PreserveNodesGPU(int nodeCount, double3 *r, double3 *r0, double3 *v, double3 *v0)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < nodeCount) {
		r0[i] = r[i];
		v0[i] = v[i];
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    ResetNodesGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void ResetNodesGPU(int nodeCount, double3 *v, double3 *v0)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < nodeCount) {
		v[i] = v0[i];
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    ResetForceGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void ResetForceGPU(int nodeCount, double3 *f)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < nodeCount) {
		f[i].x = 0.0;
		f[i].y = 0.0;
		f[i].z = 0.0;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    MobilityDragGPU
 *
 *-----------------------------------------------------------------------*/
__device__ void MobilityDragGPU(double3 dr, double L, double3 b, double &B)
{
	double3 ub;
	double  magb, dangle, Mob;
	
	if (L < 1.e-10) {
		B = 1.0; // avoid division-by-0 later
		return;
	}
	
	magb = sqrt(b.x*b.x + b.y*b.y + b.z*b.z);
	ub.x = b.x / magb;
	ub.y = b.y / magb;
	ub.z = b.z / magb;
	
	dangle = fabs(ub.x*dr.x + ub.y*dr.y + ub.z*dr.z);
	Mob = MobEdge+(MobScrew-MobEdge)*dangle;
	B = 0.5*L/Mob;
}

/*------------------------------------------------------------------------
 *
 *      Function:    SelfForceGPU
 *
 *-----------------------------------------------------------------------*/
template <unsigned int coreOnly>
__device__ void SelfForceGPU(double3 dr, double L, double3 b, double3 &f1, double3 &f2)
{
	double3 fs;
	double  Ecore, eps;

	if (coreOnly) {
		Ecore = 0.5 * TensionFactor * MU;
		eps = 1.e-06;
	} else {
		Ecore = coreEnergy;
		eps = 1.e-20;
	}
	
	if (L*L < eps) {
		return;
	}

	double tx, ty, tz, La, S;
	double bs, bs2, bex, bey, bez, be2, fL, ft;

	tx = dr.x / L;
	ty = dr.y / L;
	tz = dr.z / L;

	bs = b.x*tx + b.y*ty + b.z*tz;
	bex = b.x-bs*tx; bey = b.y-bs*ty; bez=b.z-bs*tz;
	be2 = (bex*bex+bey*bey+bez*bez);
	bs2 = bs*bs;

	La = sqrt(L*L+a*a);

	if (coreOnly) {
		S = 0.0;
	} else {
		S = (-(2*NU*La+(1-NU)*a*a/La-(1+NU)*a)/L +
			 (NU*log((La+L)/a)-(1-NU)*0.5*L/La))*MU/4/M_PI/(1-NU)*bs;
	}

	/* Ecore = MU/(4*pi) log(a/a0) */
	fL = -Ecore*(bs2+be2/(1-NU));
	ft =  Ecore*2*bs*NU/(1-NU); 

	fs.x = bex*(S+ft) + fL*tx;
	fs.y = bey*(S+ft) + fL*ty;
	fs.z = bez*(S+ft) + fL*tz;

	f2.x += fs.x;
	f2.y += fs.y;
	f2.z += fs.z;

	f1.x -= fs.x;
	f1.y -= fs.y;
	f1.z -= fs.z;
}

/*------------------------------------------------------------------------
 *
 *      Function:    ExtForceGPU
 *
 *-----------------------------------------------------------------------*/
__device__ void ExtForceGPU(double3 dr, double3 b, double *extStress, double3 &f1, double3 &f2)
{
	double3 sigb, ft;
	
	sigb.x = extStress[0]*b.x + extStress[5]*b.y + extStress[4]*b.z;
	sigb.y = extStress[5]*b.x + extStress[1]*b.y + extStress[3]*b.z;
	sigb.z = extStress[4]*b.x + extStress[3]*b.y + extStress[2]*b.z;
	
	ft.x = sigb.y*dr.z - sigb.z*dr.y;
	ft.y = sigb.z*dr.x - sigb.x*dr.z;
	ft.z = sigb.x*dr.y - sigb.y*dr.x;
	
	f1.x += 0.5*ft.x;
	f1.y += 0.5*ft.y;
	f1.z += 0.5*ft.z;
		
	f2.x += 0.5*ft.x;
	f2.y += 0.5*ft.y;
	f2.z += 0.5*ft.z;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     RemForceGPU
 *
 *-------------------------------------------------------------------------*/
__device__ void RemForceGPU(double3 dr, double3 b, double3 *sigb, double3 &f1, double3 &f2)
{
	int       i, numPoints;
	double    positions[3], weights[3];
	double    temp, mult1, mult2;
	double    sigbx, sigby, sigbz;
	double    fLinvx, fLinvy, fLinvz;
	double    pspanx, pspany, pspanz;
	
	numPoints = 3;
	positions[0] = -0.774596669241483;
	positions[1] = 0.0;
	positions[2] = -positions[0];
	weights[0] = 0.5*5.0/9.0;
	weights[1] = 0.5*8.0/9.0;
	weights[2] = weights[0];
	
	pspanx = 0.5 * dr.x;
	pspany = 0.5 * dr.y;
	pspanz = 0.5 * dr.z;
	
	for (i = 0; i < numPoints; i++) {
		
		sigbx = sigb[i].x;
		sigby = sigb[i].y;
		sigbz = sigb[i].z;

		fLinvx = (sigby*pspanz-sigbz*pspany);
		fLinvy = (sigbz*pspanx-sigbx*pspanz);
		fLinvz = (sigbx*pspany-sigby*pspanx);

		temp = weights[i]*positions[i];
		mult1 = weights[i]+temp;

		f2.x += fLinvx*mult1;
		f2.y += fLinvy*mult1;
		f2.z += fLinvz*mult1;

		mult2 = weights[i]-temp;

		f1.x += fLinvx*mult2;
		f1.y += fLinvy*mult2;
		f1.z += fLinvz*mult2;    
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    SegForceDragGPU
 *
 *-----------------------------------------------------------------------*/
template <unsigned int elasticinteraction>
__global__ void SegForceDragGPU(int segCount, double3 *r, int2 *s, double3 *b, double *extStress, double3 *fseg, double *Bseg, double3 *fmm)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < segCount) {
		
		int2 n = s[i];

		double3 r1, r2, dr, bs, f1, f2;
		double  L, B;

		r1 = r[n.x];
		r2 = r[n.y];
		bs = b[i];
		
		dr.x = r2.x - r1.x;
		dr.y = r2.y - r1.y;
		dr.z = r2.z - r1.z;
		ZImageGPU(dr);
		L = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);

		f1.x = 0.0;
		f1.y = 0.0;
		f1.z = 0.0;

		f2.x = 0.0;
		f2.y = 0.0;
		f2.z = 0.0;

		switch(elasticinteraction) {
			case 0:
				SelfForceGPU<1>(dr, L, bs, f1, f2);
				ExtForceGPU(dr, bs, extStress, f1, f2);
				break;
			case 1:
				SelfForceGPU<0>(dr, L, bs, f1, f2);
				ExtForceGPU(dr, bs, extStress, f1, f2);
				
				int numPoints = 3;
				double3 sigb[3];
				sigb[0] = fmm[i*numPoints+0];
				sigb[1] = fmm[i*numPoints+1];
				sigb[2] = fmm[i*numPoints+2];
				
				double3 fs1, fs2;
				fs1.x = 0.0; fs1.y = 0.0; fs1.z = 0.0;
				fs2.x = 0.0; fs2.y = 0.0; fs2.z = 0.0;
				
				RemForceGPU(dr, bs, sigb, fs1, fs2);
				
				f1.x += fs1.x;
				f1.y += fs1.y;
				f1.z += fs1.z;
				
				f2.x += fs2.x;
				f2.y += fs2.y;
				f2.z += fs2.z;
				break;
		}

		fseg[i*2+0] = f1;
		fseg[i*2+1] = f2;

		MobilityDragGPU(dr, L, bs, B);
		Bseg[i] = B;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    SegDragGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void SegDragGPU(int segCount, double3 *r, int2 *s, double3 *b, double *Bseg)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < segCount) {
		
		int2 n = s[i];

		double3 r1, r2, dr, bs;
		double  L, B;

		r1 = r[n.x];
		r2 = r[n.y];
		bs = b[i];
		
		dr.x = r2.x - r1.x;
		dr.y = r2.y - r1.y;
		dr.z = r2.z - r1.z;
		ZImageGPU(dr);
		L = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);

		MobilityDragGPU(dr, L, bs, B);
		Bseg[i] = B;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    SegSegForceGPU
 *
 *-----------------------------------------------------------------------*/
template <unsigned int checkFlags>
__global__ void SegSegForceGPU(int segSegCount, double3 *r, int2 *s, double3 *b, int2 *g, double3 *cc, int *gflag, double3 *fseg)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < segSegCount) {
		
		int2 segs = g[i];
		
		int2 s1 = s[segs.x];
		int2 s2 = s[segs.y];
		
		double3 b1, b2;
		b1 = b[segs.x];
		b2 = b[segs.y];
		
		double3 r1, r2, r3, r4;
		double3 f1, f2, f3, f4;
		
		int n1 = s1.x;
		int n2 = s1.y;
		int n3 = s2.x;
		int n4 = s2.y;
		
		r1 = r[n1];
		r2 = r[n2];
		r3 = r[n3];
		r4 = r[n4];
		
		int compute = 1;
		if (checkFlags == 1) {
			compute = gflag[i];
			if (compute == 0) compute = 1;
			else compute = 0;
		}
		
		double3 dr1, dr2;
		double  L1s, L2s;
		dr1.x = r2.x - r1.x;
		dr1.y = r2.y - r1.y;
		dr1.z = r2.z - r1.z;
		ZImageGPU(dr1);
		L1s = dr1.x*dr1.x + dr1.y*dr1.y + dr1.z*dr1.z;
		
		dr2.x = r4.x - r3.x;
		dr2.y = r4.y - r3.y;
		dr2.z = r4.z - r3.z;
		ZImageGPU(dr2);
		L2s = dr2.x*dr2.x + dr2.y*dr2.y + dr2.z*dr2.z;
		
		if (L1s < 1.e-20 || L2s < 1.e-20 || !compute) {
			
			f1.x = 0.0; f1.y = 0.0; f1.z = 0.0;
			f2.x = 0.0; f2.y = 0.0; f2.z = 0.0;
			f3.x = 0.0; f3.y = 0.0; f3.z = 0.0;
			f4.x = 0.0; f4.y = 0.0; f4.z = 0.0;
			
		} else {
		
			r2.x = r1.x + dr1.x;
			r2.y = r1.y + dr1.y;
			r2.z = r1.z + dr1.z;
			
			// Cell center here??
			double3 rc;
			//rc = r1;
			rc = cc[n1];
			dr1.x = r3.x - rc.x;
			dr1.y = r3.y - rc.y;
			dr1.z = r3.z - rc.z;
			ZImageGPU(dr1);
			
			r3.x = rc.x + dr1.x;
			r3.y = rc.y + dr1.y;
			r3.z = rc.z + dr1.z;
			
			dr2.x = r4.x - r3.x;
			dr2.y = r4.y - r3.y;
			dr2.z = r4.z - r3.z;
			ZImageGPU(dr2);
			
			r4.x = r3.x + dr2.x;
			r4.y = r3.y + dr2.y;
			r4.z = r3.z + dr2.z;
			
			SegSegForceIsotropicGPU(r1, r2, r3, r4, b1, b2, f1, f2, f3, f4);
			
		}
		
		fseg[i*4+0] = f1;
		fseg[i*4+1] = f2;
		fseg[i*4+2] = f3;
		fseg[i*4+3] = f4;
	}
	
	return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    ReduceForceGPU
 *
 *-----------------------------------------------------------------------*/
template <unsigned int incrementForces>
__global__ void ReduceForceGPU(int nodeCount, int2 *gpos, int *gind, double3 *fseg, double3 *f)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < nodeCount) {
		
		int2 pos = gpos[i];
		int j, k = 0;
		
		double3 fs, ft;
		ft.x = 0.0;
		ft.y = 0.0;
		ft.z = 0.0;
		
		for (k = pos.x; k < pos.y; k++) {
			j = gind[k];
			fs = fseg[j];
			ft.x += fs.x;
			ft.y += fs.y;
			ft.z += fs.z;
		}
		
		if (incrementForces == 0) {
			f[i] = ft;
		} else {
			f[i].x += ft.x;
			f[i].y += ft.y;
			f[i].z += ft.z;
		}
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    ResetArmForceGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void ResetArmForceGPU(int armCount, double3 *farms)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < armCount) {
		farms[i].x = 0.0;
		farms[i].y = 0.0;
		farms[i].z = 0.0;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    ReduceArmForceGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void ReduceArmForceGPU(int armCount, int *garms, double3 *fseg, double3 *farms)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < armCount) {
		
		int pos = garms[i];
		
		double3 fs;
		fs = fseg[pos];
		
		farms[i].x += fs.x;
		farms[i].y += fs.y;
		farms[i].z += fs.z;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    ReduceArmForceGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void ReduceArmForceGPU(int armCount, int2 *garms_pos, int *garms_ind, double3 *fseg, double3 *farms)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < armCount) {
		
		int2 pos = garms_pos[i];
		int j, k;
		
		double3 fs, ft;
		ft.x = 0.0;
		ft.y = 0.0;
		ft.z = 0.0;
		
		for (k = pos.x; k < pos.y; k++) {
			j = garms_ind[k];
			fs = fseg[j];
			ft.x += fs.x;
			ft.y += fs.y;
			ft.z += fs.z;
		}
		
		farms[i].x += ft.x;
		farms[i].y += ft.y;
		farms[i].z += ft.z;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    ReduceDragGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void ReduceDragGPU(int nodeCount, int2 *gpos, int *gind, double *Bseg, double *B)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < nodeCount) {
		
		int2 pos = gpos[i];
		int j, k = 0;
		
		double Bs, Bt;
		Bt = 0.0;
		
		for (k = pos.x; k < pos.y; k++) {
			j = gind[k]/2;
			Bs = Bseg[j];
			Bt += Bs;
		}
		B[i] = Bt;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    NodeVelocityGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void NodeVelocityGPU(int nodeCount, double3 *f, double *B, double *mob, int *n, double3 *v)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < nodeCount) {
		
		double3 fn, vn;
		double Bn, mobn[9];
		
		fn = f[i];
		Bn = B[i];
		if (Bn == 0.0) Bn = 1.0; // avoid division-by-zero
		
		for (int j = 0; j < 9; j++) {
			mobn[j] = mob[i*9+j];
		}
		
		fn.x /= Bn;
		fn.y /= Bn;
		fn.z /= Bn;
		
		// Project velocity onto glide constraints
		vn.x = mobn[0] * fn.x + mobn[1] * fn.y + mobn[2] * fn.z;
		vn.y = mobn[3] * fn.x + mobn[4] * fn.y + mobn[5] * fn.z;
		vn.z = mobn[6] * fn.x + mobn[7] * fn.y + mobn[8] * fn.z;

		// Oscillating node
		int a = n[i];
		if (a == 0) {
			vn.x = 0.0;
			vn.y = 0.0;
			vn.z = 0.0;
		}

		v[i] = vn;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    RKFStepGPU
 *
 *-----------------------------------------------------------------------*/
//template <unsigned int step>
template <unsigned int moveNodes, unsigned int blockSize>
__global__ void RKFStepGPU(int step, int nodeCount, double currDT, double3 *v, double3 *rkf, double3 *r0, double3 *r, double *e1, double *e2, double3 *f, int *n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int jmax = step+1;
	double3 vs, rrkf[6];
	double3 rold, rnew, rcur;
	
	if (i < nodeCount) {

		vs = v[i];
		rrkf[step] = vs;
		if (step < 5) rkf[i*5+step] = vs;
		for (int j = 0; j < step; j++) {
			rrkf[j] = rkf[i*5+j];
		}

		double c[6];
		if (step == 0) {
			c[0] = 1.0/4.0;
		}
		if (step == 1) {
			c[0] = 3.0/32.0;
			c[1] = 9.0/32.0;
		}
		if (step == 2) {
			c[0] = 1932.0/2197.0;
			c[1] = -7200.0/2197.0;
			c[2] = 7296.0/2197.0;
		}
		if (step == 3) {
			c[0] = 439.0/216.0;
			c[1] = -8.0;
			c[2] = 3680.0/513.0;
			c[3] = -845.0/4104.0;
		}
		if (step == 4) {
			c[0] = -8.0/27.0;
			c[1] = 2.0;
			c[2] = -3544.0/2565.0;
			c[3] = 1859.0/4104.0;
			c[4] = -11.0/40.0;
		}
		if (step == 5) {
			c[0] = 16.0/135.0;
			c[1] = 0.0;
			c[2] = 6656.0/12825.0;
			c[3] = 28561.0/56430.0;
			c[4] = -9.0/50.0;
			c[5] = 2.0/55.0;
		}

		rold = r0[i];
		rnew.x = 0.0;
		rnew.y = 0.0;
		rnew.z = 0.0;
		for (int j = 0; j < jmax; j++) {
			rnew.x += c[j]*rrkf[j].x;
			rnew.y += c[j]*rrkf[j].y;
			rnew.z += c[j]*rrkf[j].z;
		}
		rnew.x = rold.x + currDT*rnew.x;
		rnew.y = rold.y + currDT*rnew.y;
		rnew.z = rold.z + currDT*rnew.z;
		
		FoldBoxGPU(rnew);
		
		if (step == 5) {
			rcur = r[i]; // WARNING
		}
		//rcur = rnew;
		
		r[i] = rnew;
	}

	// Error calculation
	if (step == 5) {
		
		double er[6];
		er[0] =  1.0/360.0;
		er[1] =  0.0;
		er[2] = -128.0/4275.0;
		er[3] = -2197.0/75240.0;
		er[4] =  1.0/50.0;
		er[5] =  2.0/55.0;

		double errnet = 0.0;
		double relerrnet = 0.0;
		if (i < nodeCount) {
			
			double3 err;
			err.x = 0.0;
			err.y = 0.0;
			err.z = 0.0;
			for (int j = 0; j < jmax; j++) {
				err.x += er[j]*rrkf[j].x;
				err.y += er[j]*rrkf[j].y;
				err.z += er[j]*rrkf[j].z;
			}
			err.x *= currDT;
			err.y *= currDT;
			err.z *= currDT;
			errnet = sqrt(err.x*err.x+err.y*err.y+err.z*err.z); //sqrtf
			
			double3 dr;
			dr.x = rcur.x - rold.x;
			dr.y = rcur.y - rold.y;
			dr.z = rcur.z - rold.z;
			ZImageGPU(dr);
			
			double drn;
			drn = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
			if (errnet > rTolth) {
				if (drn > rTolth/rTolrel) {
					relerrnet = errnet / drn;
				} else {
					relerrnet = 2*rTolrel;
				}
			}
			
			if (moveNodes == 1) {
				if (errnet < rTol && (errnet < rTolth || errnet/drn < rTolrel)) {
					n[i] = 1; // unflag node
				} else {
					n[i] = 2; // flag node
				}
			}
		}

		// Reduce block error
		__shared__ double errtmp1[blockSize];
		__shared__ double errtmp2[blockSize];
		
		int tid = threadIdx.x;
		errtmp1[tid] = errnet;
		errtmp2[tid] = relerrnet;
		__syncthreads();

		for (unsigned int s = blockDim.x/2; s >= 1; s = s/2) {
			if (tid < s) {
				if (errtmp1[tid] < errtmp1[tid + s])
					errtmp1[tid] = errtmp1[tid + s];
				if (errtmp2[tid] < errtmp2[tid + s])
					errtmp2[tid] = errtmp2[tid + s];
			}
			__syncthreads();
		}

		if (tid == 0) {
			e1[blockIdx.x] = errtmp1[0];
			e2[blockIdx.x] = errtmp2[0];
		}
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    FlagNodesGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void FlagNodesGPU(int nodeCount, int *n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < nodeCount) {
		n[i] = 1;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    MoveInteractionsGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void MoveInteractionsGPU(int segSegCount, double rgs, int2 *s, int2 *g, double *gdist2, int *n, int *gflag)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < segSegCount) {
		
		int2 segs = g[i];
		int2 s1 = s[segs.x];
		int2 s2 = s[segs.y];
		
		int n1 = s1.x;
		int n2 = s1.y;
		int n3 = s2.x;
		int n4 = s2.y;
		
		int a1 = n[n1];
		int a2 = n[n2];
		int a3 = n[n3];
		int a4 = n[n4];
		
		// Flag the interaction to be moved to group 4 
		// if any of its nodes is flagged 2
		if ((a1-2)*(a2-2)*(a3-2)*(a4-2) == 0) {
			if (gdist2[i] <= rgs && gflag[i] == 0) {
				gflag[i] = 4;
			}
		}
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    ForwardProgressCheckGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void ForwardProgressCheckGPU(int nodeCount, double3 *v0, double3 *v, int *n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < nodeCount) {
		double3 vold, vnew;
		vold = v0[i];
		vnew = v[i];
		
		double vv;
		vv = vold.x*vnew.x + vold.y*vnew.y + vold.z*vnew.z;
		if (vv < 0.0) n[i] = 0;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    RKFIntegratorGPU
 *
 *-----------------------------------------------------------------------*/
void RKFIntegratorGPU(Home_t *home, Device_t *device, int reqType, int nSubcyc)
{
		int     i, threadsPerBlock;
		int     blocksNodes, blocksSegs, blocksSegSegs, blocksArms; 
		int     nodeCount, segCount, segSegCount, armCount;
		double  errMax, relErrMax, newDT, rg9s, tmp_rg9s;
		int2    *seg;
		int2    *group;
		int2    *spos, *gpos;
		int     *gind, *sind, *gflag;
		double3 *fseg;
		double3 *fmm;
		Param_t *param;
		
		param = home->param;
		
		
		/* Segment list */
		seg  = device->g1;
		spos = device->g1pos;
		sind = device->g1ind;
		segCount   = device->segCount;
		blocksSegs = device->blocksSegs;
		armCount   = device->armCount;
		blocksArms = device->blocksArms;
		threadsPerBlock = device->threadsPerBlock;
		
		fmm = NULL;
		gflag = NULL;
		
		if (reqType == FULL) {
			
			nodeCount   = device->nodeCount;
			blocksNodes = device->blocksNodes;
			
			fseg = device->fseg;
			gpos = device->g1pos;
			gind = device->g1ind;
			
			if (param->elasticinteraction) {
				fmm  = device->fmm;
			}
			
			if (nSubcyc == -1) {
				gflag = device->g0flag;
			}
			
		} else if (reqType == GROUP0) {
			
			nodeCount     = device->nodeCount;
			blocksNodes   = device->blocksNodes;
			segSegCount   = device->nSegSeg0;
			blocksSegSegs = device->blocksSegSegs0;
			
			fseg  = device->f0;
			group = device->g0;
			gpos  = device->g0pos;
			gind  = device->g0ind;
			gflag = device->g0flag;
			
			rg9s = MAX(MAX(MAX(param->rg1,param->rg2),param->rg3),param->rg4) * 2;
			rg9s = rg9s * rg9s;
			rg9s = 5000.0 *5000.0 ; // WARNING
			//This should probably not be set to a constant ????
			
		} else if (reqType == GROUP1) {
			
			nodeCount   = device->nodeCount;
			blocksNodes = device->blocksNodes;
			
			fseg = device->fseg;
			gpos = device->g1pos;
			gind = device->g1ind;
			
			if (param->forceCutOff == 0 && param->fmEnabled) {
				fmm  = device->fmm;
			}
			
		} else if (reqType == GROUP2) {
			
			nodeCount     = device->nodeCount;
			blocksNodes   = device->blocksNodes;
			segSegCount   = device->nSegSeg2;
			blocksSegSegs = device->blocksSegSegs2;
			
			fseg  = device->f2;
			group = device->g2;
			gpos  = device->g2pos;
			gind  = device->g2ind;
			
		} else if (reqType == GROUP3) {
			
			nodeCount     = device->nodeCount;
			blocksNodes   = device->blocksNodes;
			segSegCount   = device->nSegSeg3;
			blocksSegSegs = device->blocksSegSegs3;
			
			fseg  = device->f3;
			group = device->g3;
			gpos  = device->g3pos;
			gind  = device->g3ind;
			
		} else if (reqType == GROUP4) {
			
			nodeCount     = device->nodeCount;
			blocksNodes   = device->blocksNodes;
			segSegCount   = device->nSegSeg4;
			blocksSegSegs = device->blocksSegSegs4;
			
			fseg  = device->f4;
			group = device->g4;
			gpos  = device->g4pos;
			gind  = device->g4ind;
			
		} else {
			Fatal("GPU subcycling is not available for this subGroup yet!");
		}
		
		
		/* Flag all nodes for subcycling */
		if (nSubcyc <= 0) {
			FlagNodesGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, device->n);
		}
		
		switch (param->elasticinteraction) {
			case 0:
				if (reqType == FULL || reqType == GROUP1) {
					SegForceDragGPU<0><<<blocksSegs,threadsPerBlock>>>(segCount, device->r, seg, device->b, device->esig, fseg, device->Bseg, fmm);
				} else {
					/* Dummy, just to avoid unitialized memory read later */
					SegDragGPU<<<blocksSegs,threadsPerBlock>>>(segCount, device->r, seg, device->b, device->Bseg);
				}
				break;
			case 1:
				if (reqType == FULL || reqType == GROUP1) {
					SegForceDragGPU<1><<<blocksSegs,threadsPerBlock>>>(segCount, device->r, seg, device->b, device->esig, fseg, device->Bseg, fmm);
					if (reqType == FULL) {
						if (nSubcyc == -1) {
							SegSegForceGPU<1><<<device->blocksSegSegs0,threadsPerBlock>>>(device->nSegSeg0, device->r, seg, device->b, device->g0, device->cc, gflag, device->f0);
						} else {
							SegSegForceGPU<0><<<device->blocksSegSegs0,threadsPerBlock>>>(device->nSegSeg0, device->r, seg, device->b, device->g0, device->cc, gflag, device->f0);
						}
						SegSegForceGPU<0><<<device->blocksSegSegs2,threadsPerBlock>>>(device->nSegSeg2, device->r, seg, device->b, device->g2, device->cc, gflag, device->f2);
						SegSegForceGPU<0><<<device->blocksSegSegs3,threadsPerBlock>>>(device->nSegSeg3, device->r, seg, device->b, device->g3, device->cc, gflag, device->f3);
						SegSegForceGPU<0><<<device->blocksSegSegs4,threadsPerBlock>>>(device->nSegSeg4, device->r, seg, device->b, device->g4, device->cc, gflag, device->f4);
					}
				} else {
					if (reqType == GROUP0) {
						SegSegForceGPU<1><<<blocksSegSegs,threadsPerBlock>>>(segSegCount, device->r, seg, device->b, group, device->cc, gflag, fseg);
					} else {
						SegSegForceGPU<0><<<blocksSegSegs,threadsPerBlock>>>(segSegCount, device->r, seg, device->b, group, device->cc, gflag, fseg);
					}
					SegDragGPU<<<blocksSegs,threadsPerBlock>>>(segCount, device->r, seg, device->b, device->Bseg);
				}
				break;
		}
		
		if (reqType == FULL) {
			ResetForceGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, device->f);
			ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g0pos, device->g0ind, device->f0, device->f);
			ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g1pos, device->g1ind, device->fseg, device->f);
			ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g2pos, device->g2ind, device->f2, device->f);
			ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g3pos, device->g3ind, device->f3, device->f);
			ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g4pos, device->g4ind, device->f4, device->f);
		} else {
			ReduceForceGPU<0><<<blocksNodes,threadsPerBlock>>>(nodeCount, gpos, gind, fseg, device->f);
		}
		ReduceDragGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, spos, sind, device->Bseg, device->B);
		
		NodeVelocityGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, device->f, device->B, device->mob, device->n, device->v);
		PreserveNodesGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, device->r, device->r0, device->v, device->v0);
		
		
		/* Force velocity calculation only */
		if (nSubcyc == -1) {
			if (reqType == FULL) {
				/* When we are doing a full calculation, we need to communicate the arm forces back as well */
				ResetArmForceGPU<<<blocksArms,threadsPerBlock>>>(armCount, device->farms);
				ReduceArmForceGPU<<<blocksArms,threadsPerBlock>>>(armCount, device->g0arms_pos, device->g0arms_ind, device->f0, device->farms);
				ReduceArmForceGPU<<<blocksArms,threadsPerBlock>>>(armCount, device->g1arms, device->fseg, device->farms);
				ReduceArmForceGPU<<<blocksArms,threadsPerBlock>>>(armCount, device->g2arms_pos, device->g2arms_ind, device->f2, device->farms);
				ReduceArmForceGPU<<<blocksArms,threadsPerBlock>>>(armCount, device->g3arms_pos, device->g3arms_ind, device->f3, device->farms);
				ReduceArmForceGPU<<<blocksArms,threadsPerBlock>>>(armCount, device->g4arms_pos, device->g4arms_ind, device->f4, device->farms);
			}
			return;
		}
		
		
		/* Grab the time step */
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
		
		/* RKF integration: initialize convergence loop */
		int convergent =  0;
		int incrDelta  =  1;
		int iTry       = -1;
		
		/* If there are no interactions, just skip this loop */
		if (reqType == GROUP0 && segSegCount == 0) {
			convergent = 1;
			errMax = 0.0;
		}

		while (!convergent) {
			iTry++;

			for (i = 0; i < 5; i++) {
				
				RKFStepGPU<0,1><<<blocksNodes,threadsPerBlock>>>(i, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n);
				
				switch (param->elasticinteraction) {
					case 0:
						if (reqType == FULL || reqType == GROUP1) {
							SegForceDragGPU<0><<<blocksSegs,threadsPerBlock>>>(segCount, device->r, seg, device->b, device->esig, fseg, device->Bseg, fmm);
						}
						break;
					case 1:
						if (reqType == FULL || reqType == GROUP1) {
							SegForceDragGPU<1><<<blocksSegs,threadsPerBlock>>>(segCount, device->r, seg, device->b, device->esig, fseg, device->Bseg, fmm);
							if (reqType == FULL) {
								SegSegForceGPU<0><<<device->blocksSegSegs0,threadsPerBlock>>>(device->nSegSeg0, device->r, seg, device->b, device->g0, device->cc, gflag, device->f0);
								SegSegForceGPU<0><<<device->blocksSegSegs2,threadsPerBlock>>>(device->nSegSeg2, device->r, seg, device->b, device->g2, device->cc, gflag, device->f2);
								SegSegForceGPU<0><<<device->blocksSegSegs3,threadsPerBlock>>>(device->nSegSeg3, device->r, seg, device->b, device->g3, device->cc, gflag, device->f3);
								SegSegForceGPU<0><<<device->blocksSegSegs4,threadsPerBlock>>>(device->nSegSeg4, device->r, seg, device->b, device->g4, device->cc, gflag, device->f4);
							}
						} else {
							if (reqType == GROUP0) {
								SegSegForceGPU<1><<<blocksSegSegs,threadsPerBlock>>>(segSegCount, device->r, seg, device->b, group, device->cc, gflag, fseg);
							} else {
								SegSegForceGPU<0><<<blocksSegSegs,threadsPerBlock>>>(segSegCount, device->r, seg, device->b, group, device->cc, gflag, fseg);
							}
							SegDragGPU<<<blocksSegs,threadsPerBlock>>>(segCount, device->r, seg, device->b, device->Bseg);
						}
						break;
				}
				
				if (reqType == FULL) {
					ResetForceGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, device->f);
					ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g0pos, device->g0ind, device->f0, device->f);
					ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g1pos, device->g1ind, device->fseg, device->f);
					ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g2pos, device->g2ind, device->f2, device->f);
					ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g3pos, device->g3ind, device->f3, device->f);
					ReduceForceGPU<1><<<blocksNodes,threadsPerBlock>>>(nodeCount, device->g4pos, device->g4ind, device->f4, device->f);
				} else {
					ReduceForceGPU<0><<<blocksNodes,threadsPerBlock>>>(nodeCount, gpos, gind, fseg, device->f);
				}
				ReduceDragGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, spos, sind, device->Bseg, device->B);
				
				NodeVelocityGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, device->f, device->B, device->mob, device->n, device->v);
			}
			
			if (reqType == GROUP0 && iTry < param->nTry) {
				switch (threadsPerBlock) {
					case 1024: RKFStepGPU<1,1024><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					case 512:  RKFStepGPU<1,512><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					case 256:  RKFStepGPU<1,256><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					case 128:  RKFStepGPU<1,128><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					case 64:   RKFStepGPU<1,64><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					default:   Fatal("RKFStepGPU is not defined for threadsPerBlock = %d", threadsPerBlock); break;
				}
			} else {
				switch (threadsPerBlock) {
					case 1024: RKFStepGPU<0,1024><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					case 512:  RKFStepGPU<0,512><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					case 256:  RKFStepGPU<0,256><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					case 128:  RKFStepGPU<0,128><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					case 64:   RKFStepGPU<0,64><<<blocksNodes,threadsPerBlock>>>(5, nodeCount, newDT, device->v, device->rkf, device->r0, device->r, device->e1, device->e2, device->f, device->n); break;
					default:   Fatal("RKFStepGPU is not defined for threadsPerBlock = %d", threadsPerBlock); break;
				}
			}
			
			/* Calculate maximum error */
#if 1
			errMax = *(thrust::max_element(device->e1_ptr, device->e1_ptr + blocksNodes));
			relErrMax = *(thrust::max_element(device->e2_ptr, device->e2_ptr + blocksNodes));
#else
			double *e1 = (double*)malloc(sizeof(double)*blocksNodes); //ALLOCATE this once per time-step if used
			double *e2 = (double*)malloc(sizeof(double)*blocksNodes);
			HANDLE_ERROR(cudaMemcpy(e1, device->e1, sizeof(double)*blocksNodes, cudaMemcpyDeviceToHost));
			HANDLE_ERROR(cudaMemcpy(e2, device->e2, sizeof(double)*blocksNodes, cudaMemcpyDeviceToHost));
			errMax = 0.0;
			relErrMax = 0.0;
			for (i = 0; i < blocksNodes; i++) {
				errMax = MAX(e1[i], errMax);
				relErrMax = MAX(e2[i], relErrMax);
			}
			free(e1);
			free(e2);
#endif
			
			if (errMax < param->rTol && relErrMax < param->rTolrel) {
				convergent = 1;
				
				/* Flag oscillating nodes for subsequent cycles */
				if (reqType > GROUP0 && nSubcyc > 3) {
					ForwardProgressCheckGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, device->v0, device->v, device->n);
				}
				
			} else {
				
				if ((iTry < param->nTry) && (reqType == GROUP0)) {
					
					tmp_rg9s = rg9s * (iTry+1) * (iTry+1);
					MoveInteractionsGPU<<<blocksSegSegs,threadsPerBlock>>>(segSegCount, tmp_rg9s, seg, group, device->g0dist2, device->n, gflag);
				
				} else {
				
					// Restart with old velocities
					ResetNodesGPU<<<blocksNodes,threadsPerBlock>>>(nodeCount, device->v, device->v0);
					
					incrDelta = 0;
					newDT    *= param->dtDecrementFact;

					if ((newDT < 1.0e-20) && (home->myDomain == 0)) {
						Fatal("RKFIntegrator(): Timestep has dropped below\n"
							  "minimal threshold to %e.  Aborting!", newDT);
					}
				
				}
			}

		}
		
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
				tmp2 = errMax/param->rTol;
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
		
}

/*------------------------------------------------------------------------
 *
 *      Function:    GetMinDistGPU
 *      Description: Calculate the minimum distance between two segments
 *
 *-----------------------------------------------------------------------*/
__device__ void GetMinDistGPU(double3 r1, double3 r2, double3 r3, double3 r4, double &dist2)
{
	int     i, pos;
	int     icase, didDist2;
	double  A, B, C, D, E;
	double  eps = 1.0e-12;
	double  distx, disty, distz, d2, d2min;
	double  dist[4], L1, L2;

	double3 seg1L, seg2L;
	double3 r1mr3, r2mr1, r4mr3, r4mr1, r3mr2, r4mr2;
	double  M[2][2], rhs[2], sol[2], detM;
	double  trial[4][2];
	
	r1mr3.x = r1.x - r3.x;
	r1mr3.y = r1.y - r3.y;
	r1mr3.z = r1.z - r3.z;
	
	r2mr1.x = r2.x - r1.x;
	r2mr1.y = r2.y - r1.y;
	r2mr1.z = r2.z - r1.z;
	
	r4mr3.x = r4.x - r3.x;
	r4mr3.y = r4.y - r3.y;
	r4mr3.z = r4.z - r3.z;
	
	seg1L = r2mr1;
	seg2L = r4mr3;
	
	M[0][0] = DotProductGPU(r2mr1, r2mr1);
	M[1][0] =-DotProductGPU(r4mr3, r2mr1);
	M[1][1] = DotProductGPU(r4mr3, r4mr3);
	M[0][1] = M[1][0];
	
	rhs[0] = -DotProductGPU(r2mr1, r1mr3);
	rhs[1] =  DotProductGPU(r4mr3, r1mr3);
	
	detM = 1.0 - M[1][0] * M[1][0] / M[0][0] / M[1][1];

	A = M[0][0];
	B = -2.0 * rhs[0];
	C = -2.0 * M[1][0];
	D = -2.0 * rhs[1];
	E = M[1][1];

	didDist2 = 0;
	
/*
 *      If segment 1 is just a point...
 */
        if (A < eps) {
            L1 = 0.0;
            if (E < eps) L2 = 0.0;
            else L2 = -0.5 * D / E;

/*
 *      If segment 2 is just a point...
 */
        } else if (E < eps) {
            L2 = 0.0;
            if (A < eps) L1 = 0.0;
            else L1 = -0.5 * B / A;
/*
 *      If segments are parallel
 */
		} else if (detM<1e-6) {
			
			r4mr1.x = r4.x - r1.x;
			r4mr1.y = r4.y - r1.y;
			r4mr1.z = r4.z - r1.z;
			
			r3mr2.x = r3.x - r2.x;
			r3mr2.y = r3.y - r2.y;
			r3mr2.z = r3.z - r2.z;
			
			r4mr2.x = r4.x - r2.x;
			r4mr2.y = r4.y - r2.y;
			r4mr2.z = r4.z - r2.z;
			
			dist[0] = DotProductGPU(r1mr3, r1mr3);
            dist[1] = DotProductGPU(r4mr1, r4mr1);
            dist[2] = DotProductGPU(r3mr2, r3mr2);
            dist[3] = DotProductGPU(r4mr2, r4mr2);

            dist2 = dist[0];
            pos = 1;

            for (i = 1; i < 4; i++) {
                if (dist[i] < dist2) {
                    dist2 = dist[i];
                    pos = i+1;
                }
            }

            L1 = floor((double)pos/2.1);
			L2 = (double)(1 - (pos % 2));
			didDist2 = 1;
/*
 *		Solve the general case
 */
		} else { 
			detM *= M[0][0]*M[1][1];
			sol[0] = ( M[1][1]*rhs[0] - M[0][1]*rhs[1]) / detM;
			sol[1] = (-M[1][0]*rhs[0] + M[0][0]*rhs[1]) / detM;

			if ((sol[0]>=0) && (sol[0]<=1) && (sol[1]>=0) && (sol[1]<=1)) {
				/* we are done here */
				L1 = sol[0];
				L2 = sol[1];

			} else {

				/* enumerate four cases */
				/* alpha = 0 */
				icase = 0;
				trial[icase][0] = 0;
				trial[icase][1] = (rhs[1] - M[1][0]*trial[icase][0]) / M[1][1];

				/* alpha = 1 */
				icase = 1;
				trial[icase][0] = 1;
				trial[icase][1] = (rhs[1] - M[1][0]*trial[icase][0]) / M[1][1];

				/* beta = 0 */
				icase = 2;
				trial[icase][1] = 0;
				trial[icase][0] = (rhs[0] - M[0][1]*trial[icase][1]) / M[0][0];

				/* beta = 1 */
				icase = 3;
				trial[icase][1] = 1;
				trial[icase][0] = (rhs[0] - M[0][1]*trial[icase][1]) / M[0][0];

				/* find the minimum out of four trials */
				d2min = 1e100;
				for(icase = 0; icase < 4; icase++) {
					trial[icase][0] = min(max(trial[icase][0], 0.0), 1.0);
					trial[icase][1] = min(max(trial[icase][1], 0.0), 1.0);
					distx = r1.x + (seg1L.x * trial[icase][0]) 
					      - r3.x - (seg2L.x * trial[icase][1]);
					disty = r1.y + (seg1L.y * trial[icase][0]) 
					      - r3.y - (seg2L.y * trial[icase][1]);
					distz = r1.z + (seg1L.z * trial[icase][0]) 
					      - r3.z - (seg2L.z * trial[icase][1]);

					d2 = distx*distx + disty*disty + distz*distz;
					if (d2<d2min) {
						L1 = trial[icase][0];
						L2 = trial[icase][1];
						d2min = d2;
					}
				}
				dist2 = d2min;
				didDist2 = 1;
			}
		} 

/*
 *      Make sure L1 and L2 are between 0 and 1
 */
        L1 = min(max(L1, 0.0), 1.0);
        L2 = min(max(L2, 0.0), 1.0);

		if (!didDist2) {
			distx = r1.x + (seg1L.x * L1) - r3.x - (seg2L.x * L2);
			disty = r1.y + (seg1L.y * L1) - r3.y - (seg2L.y * L2);
			distz = r1.z + (seg1L.z * L1) - r3.z - (seg2L.z * L2);

			dist2 = distx*distx + disty*disty + distz*distz;
		}
	
	return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    GetMinDistSegSegGPU
 *      Description: Determine segment pairs interaction group based
 *                   on the segment / segment distance
 *
 *-----------------------------------------------------------------------*/
__global__ void GetMinDistSegSegGPU(int segSegCount, double3 *r, int2 *s, int2 *g, int *rg, double *gdist2, int *gflag)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < segSegCount) {
		
		double dist2;
		
		int2 segs = g[i];
		
		int2 s1 = s[segs.x];
		int2 s2 = s[segs.y];
		
		double3 r1, r2, r3, r4;
		
		int n1 = s1.x;
		int n2 = s1.y;
		int n3 = s2.x;
		int n4 = s2.y;
		
		r1 = r[n1];
		r2 = r[n2];
		r3 = r[n3];
		r4 = r[n4];
		
		double3 dr1, dr2;
		double  L1s, L2s;
		dr1.x = r2.x - r1.x;
		dr1.y = r2.y - r1.y;
		dr1.z = r2.z - r1.z;
		ZImageGPU(dr1);
		L1s = dr1.x*dr1.x + dr1.y*dr1.y + dr1.z*dr1.z;
		
		dr2.x = r4.x - r3.x;
		dr2.y = r4.y - r3.y;
		dr2.z = r4.z - r3.z;
		ZImageGPU(dr2);
		L2s = dr2.x*dr2.x + dr2.y*dr2.y + dr2.z*dr2.z;
		
		if (L1s < 1.e-20 || L2s < 1.e-20) {
			dist2 = -1.0;
		} else {
			
			r2.x = r1.x + dr1.x;
			r2.y = r1.y + dr1.y;
			r2.z = r1.z + dr1.z;
			
			dr1.x = r3.x - r1.x;
			dr1.y = r3.y - r1.y;
			dr1.z = r3.z - r1.z;
			ZImageGPU(dr1);
			
			r3.x = r1.x + dr1.x;
			r3.y = r1.y + dr1.y;
			r3.z = r1.z + dr1.z;
			
			dr2.x = r4.x - r3.x;
			dr2.y = r4.y - r3.y;
			dr2.z = r4.z - r3.z;
			ZImageGPU(dr2);
			
			r4.x = r3.x + dr2.x;
			r4.y = r3.y + dr2.y;
			r4.z = r3.z + dr2.z;
			
			int hinge = 0;
			if (n1 == n3) {
				hinge = 1;
			} else if (n2 == n3) {
				hinge = 2;
			} else if (n2 == n4) {
				hinge = 3;
			} else if (n1 == n4) {
				hinge = 4;
			}
			
			if (!hinge) {
				GetMinDistGPU(r1, r2, r3, r4, dist2);
			} else {
				
				double3 m, h1, h2;
				
				if (hinge == 1) {
					m = r1;
					h1 = r2;
					h2 = r4;
					
				} else if (hinge == 2) {
					m = r2;
					h1 = r1;
					h2 = r4;
					
				} else if (hinge == 3) {
					m = r2;
					h1 = r1;
					h2 = r3;
					
				} else if (hinge == 4) {
					m = r1;
					h1 = r2;
					h2 = r3;
				}
				
				if (L1s>L2s) {
					GetMinDistGPU(m , h1, h2, h2, dist2); 
				} else {
					GetMinDistGPU(m, h2, h1, h1, dist2); 
				}
			}
			
		}
		
		/*
		if ((*dist2 > param->cutoff2 * param->cutoff2) &&
		    (param->forceCutOff)) *dist2 = -1.0;
		*/
		
		int SubGroup;
		if (dist2 < 0) {
			SubGroup = -1;
		} else if (dist2 < rg[0]) {
			SubGroup = 1;
		} else if (dist2 < rg[1]) {
			SubGroup = 2;
		} else if (dist2 < rg[2]) {
			SubGroup = 3;
		} else if (dist2 < rg[3]) {
			SubGroup = 4;
		} else {
			SubGroup = 0;
		}
		
		gflag[i] = SubGroup;
		gdist2[i] = dist2;
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    SegSegListGPU
 *      Description: Assign each segment pair interaction to a subcycle
 *                   group based on their interaction distance
 *
 *-----------------------------------------------------------------------*/
void SegSegListGPU(Home_t *home, Device_t *device)
{
		int segSegCount, blocksSegSegs, threadsPerBlock;
		Param_t *param = home->param;
		
		segSegCount     = device->nSegSeg0;
		blocksSegSegs   = device->blocksSegSegs0;
		threadsPerBlock = device->threadsPerBlock;
		
		int2 *seg = device->g1;
		
		int rg_host[4], *rg_device;
		rg_host[0] = param->rg1 * param->rg1;
		rg_host[1] = param->rg2 * param->rg2;
		rg_host[2] = param->rg3 * param->rg3;
		rg_host[3] = param->rg4 * param->rg4;
		
		HANDLE_ERROR(cudaMalloc(&rg_device, sizeof(int)*4));
		HANDLE_ERROR(cudaMemcpy(rg_device, rg_host, sizeof(int)*4, cudaMemcpyHostToDevice));
		
		GetMinDistSegSegGPU<<<blocksSegSegs,threadsPerBlock>>>(segSegCount, device->r, seg, device->g0, rg_device, device->g0dist2, device->g0flag);
		
		HANDLE_ERROR(cudaFree(rg_device));
		
		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    PackReductionGPU
 * 		Description: Pack the reduction array used to assemble nodal
 * 					 forces from GPU forces array.
 *
 *-----------------------------------------------------------------------*/
__global__ void PackReductionGPU(int segSegCount, int2 *g, int2 *s, int *nind, int *aind, int2 *gpos, int2 *garms_pos, int *gind, int *garms_ind)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < segSegCount) {
		
		int i1, i2, i3, i4;
		int pos1, pos2, pos3, pos4;
		
		int2 segs = g[i];
		int2 s1 = s[segs.x];
		int2 s2 = s[segs.y];
		
		i1 = s1.x;
		i2 = s1.y;
		i3 = s2.x;
		i4 = s2.y;
		
		pos1 = gpos[i1].x + nind[4*i+0];
		pos2 = gpos[i2].x + nind[4*i+1];
		pos3 = gpos[i3].x + nind[4*i+2];
		pos4 = gpos[i4].x + nind[4*i+3];
		
		gind[pos1] = i*4+0;
		gind[pos2] = i*4+1;
		gind[pos3] = i*4+2;
		gind[pos4] = i*4+3;
		
		i1 = 2*segs.x+0;
		i2 = 2*segs.x+1;
		i3 = 2*segs.y+0;
		i4 = 2*segs.y+1;
		
		pos1 = garms_pos[i1].x + aind[4*i+0];
		pos2 = garms_pos[i2].x + aind[4*i+1];
		pos3 = garms_pos[i3].x + aind[4*i+2];
		pos4 = garms_pos[i4].x + aind[4*i+3];
		
		garms_ind[pos1] = i*4+0;
		garms_ind[pos2] = i*4+1;
		garms_ind[pos3] = i*4+2;
		garms_ind[pos4] = i*4+3;
		
	}
}

/*------------------------------------------------------------------------
 *
 *      Function:    SendSubgroupGPU
 *      Description: Pack, allocate and send segments / segments 
 *                   information for the different groups on the GPU.
 *
 *-----------------------------------------------------------------------*/
void SendSubgroupGPU(Home_t *home, Device_t *device, int subGroup)
{
		int         i, j, k, nodeCount, armCount, nSegSeg;
		int         s1, s2, n1, n2, n3, n4, a1, a2, a3, a4;
		Node_t      *node1, *node2, *node3, *node4;
		SegSeg_t    *SegSegList;
		Subcyc_t    *subcyc;
		
		subcyc  = home->subcyc;
		nodeCount = device->nodeCount;
		armCount = device->armCount;
		
/*
 *		Group1 (segments) packing for GPU has already be done at this stage
 */		
		if (subGroup == GROUP1) {
			return;
		}
		
/*
 *		Pack groups and segment / segment interactions for GPU
 */		
		//if (param->elasticinteraction) {
				
			if (subGroup == GROUP0) {
				SegSegList = subcyc->SegSegListG0;
				nSegSeg    = subcyc->SegSegListG0_cnt;
			} else if (subGroup == GROUP2) {
				SegSegList = subcyc->SegSegListG2;
				nSegSeg    = subcyc->SegSegListG2_cnt;
			} else if (subGroup == GROUP3) {
				SegSegList = subcyc->SegSegListG3;
				nSegSeg    = subcyc->SegSegListG3_cnt;
			} else if (subGroup == GROUP4) {
				SegSegList = subcyc->SegSegListG4;
				nSegSeg    = subcyc->SegSegListG4_cnt;
			}
			
			int cntSegSeg = 0;
			for (j = 0; j < nSegSeg; j++) {
				if (SegSegList[j].flag == 0) continue;
				cntSegSeg++;
			}
			
			// Check block size
			int blocksSegSegs = (cntSegSeg + device->threadsPerBlock - 1) / device->threadsPerBlock;
			if (blocksSegSegs > home->deviceProp->maxBlocks) {
				Fatal("Max number of CUDA blocks exceeded for group %d!", subGroup-GROUP0);
			}
			
			
			int2 *gList = (int2*)malloc(cntSegSeg*sizeof(int2));
			
			int2 *redg_pos = (int2*)malloc(nodeCount*sizeof(int2));
			int  *redg_ind = (int*)malloc(4*cntSegSeg*sizeof(int));
			
			int2 *garm_pos = (int2*)malloc(armCount*sizeof(int2));
			int  *garm_ind = (int*)malloc(4*cntSegSeg*sizeof(int));
			
			
			if (subGroup == GROUP0) {			
/*				
 * 				Initially, all segments / segments interactions are in group 0
 * 				when using the GPU subcycle integrator (this is to maximize
 * 				code performance). Pack the interactions information and store
 * 				the positions of the interactions forces to be retreived from the
 * 				GPU force array for each node in the reduction array. In the case 
 * 				of group 0, the reduction array is built on the GPU to gain time.
 */
				int npos = 0;
				int apos = 0;
				int nCount = 0;
				
				for (i = 0; i < home->newNodeKeyPtr; i++) {
					if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
					redg_pos[nCount].x = npos;
					npos += node1->numInt;
					redg_pos[nCount].y = npos;
					nCount++;
				}
				for (j = 0; j < subcyc->SegListG1_cnt; j++) {
					node1 = subcyc->SegListG1[j].seg->node1;
					node2 = subcyc->SegListG1[j].seg->node2;
					a1 = subcyc->SegListG1[j].seg->armID12;
					a2 = subcyc->SegListG1[j].seg->armID21;
					
					garm_pos[2*j+0].x = apos;
					apos += node1->armInt[a1];
					garm_pos[2*j+0].y = apos;
					node1->armInt[a1] = 0;
					
					garm_pos[2*j+1].x = apos;
					apos += node2->armInt[a2];
					garm_pos[2*j+1].y = apos;
					node2->armInt[a2] = 0;
				}
				
				cntSegSeg = 0;
				for (j = 0; j < nSegSeg; j++) {
					if (SegSegList[j].flag == 0) continue;
					gList[cntSegSeg].x = SegSegList[j].seg1->subindex;
					gList[cntSegSeg].y = SegSegList[j].seg2->subindex;
					cntSegSeg++;
				}
				
				// Allocate memory on device
				device->nSegSeg0 = cntSegSeg;
				device->blocksSegSegs0 = blocksSegSegs;
				
				HANDLE_ERROR(cudaMalloc(&device->f0, sizeof(double3)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g0, sizeof(int2)*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g0pos, sizeof(int2)*nodeCount));
				HANDLE_ERROR(cudaMalloc(&device->g0ind, sizeof(int)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g0arms_pos, sizeof(int2)*armCount));
				HANDLE_ERROR(cudaMalloc(&device->g0arms_ind, sizeof(int)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g0flag, sizeof(int)*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g0dist2, sizeof(double)*cntSegSeg));
				
				int *device_nind, *device_aind;
				HANDLE_ERROR(cudaMalloc(&device_nind, sizeof(int)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device_aind, sizeof(int)*4*cntSegSeg));
				
				HANDLE_ERROR(cudaMemcpy(device->g0, gList, sizeof(int2)*cntSegSeg, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g0pos, redg_pos, sizeof(int2)*nodeCount, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g0arms_pos, garm_pos, sizeof(int2)*armCount, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device_nind, subcyc->SegSegListNodeInd, sizeof(int)*4*cntSegSeg, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device_aind, subcyc->SegSegListArmInd, sizeof(int)*4*cntSegSeg, cudaMemcpyHostToDevice));
				
				// Reduction kernel
				PackReductionGPU<<<blocksSegSegs,device->threadsPerBlock>>>(cntSegSeg, device->g0, device->g1, device_nind, device_aind, device->g0pos, device->g0arms_pos, device->g0ind, device->g0arms_ind);
				
				HANDLE_ERROR(cudaFree(device_nind));
				HANDLE_ERROR(cudaFree(device_aind));
				
			} else {
/*			
 * 				For other groups, pack both the segment pairs and 
 * 				the reduction array on the CPU.
 */
				std::vector<std::vector<int> > redg(nodeCount);
				std::vector<std::vector<int> > redgarms(armCount);

				cntSegSeg = 0;
				for (j = 0; j < nSegSeg; j++) {
					if (SegSegList[j].flag == 0) continue;
					
					s1 = SegSegList[j].seg1->subindex;
					s2 = SegSegList[j].seg2->subindex;
					
					gList[cntSegSeg].x = s1;
					gList[cntSegSeg].y = s2;
					
					node1 = SegSegList[j].seg1->node1;
					node2 = SegSegList[j].seg1->node2;
					node3 = SegSegList[j].seg2->node1;
					node4 = SegSegList[j].seg2->node2;
					
					n1 = node1->subindex;
					n2 = node2->subindex;
					n3 = node3->subindex;
					n4 = node4->subindex;
					
					redg[n1].push_back(cntSegSeg*4+0);
					redg[n2].push_back(cntSegSeg*4+1);
					redg[n3].push_back(cntSegSeg*4+2);
					redg[n4].push_back(cntSegSeg*4+3);
					
					a1 = SegSegList[j].seg1->armID12;
					a2 = SegSegList[j].seg1->armID21;
					a3 = SegSegList[j].seg2->armID12;
					a4 = SegSegList[j].seg2->armID21;
					
					redgarms[node1->armid[a1]].push_back(cntSegSeg*4+0);
					redgarms[node2->armid[a2]].push_back(cntSegSeg*4+1);
					redgarms[node3->armid[a3]].push_back(cntSegSeg*4+2);
					redgarms[node4->armid[a4]].push_back(cntSegSeg*4+3);
					
					cntSegSeg++;
				}
			
				// Pack reduction array
				int ind = 0;
				for (j = 0; j < nodeCount; j++) {
					redg_pos[j].x = ind;
					for (k = 0; k < redg[j].size(); k++) {
						redg_ind[ind++] = redg[j][k];
					}
					redg_pos[j].y = ind;
				}
				if (ind != 4*cntSegSeg) Fatal("Group%d reduction array size error!", subGroup-GROUP0);
			
				// Pack arm reduction array
				ind = 0;
				for (j = 0; j < armCount; j++) {
					garm_pos[j].x = ind;
					for (k = 0; k < redgarms[j].size(); k++) {
						garm_ind[ind++] = redgarms[j][k];
					}
					garm_pos[j].y = ind;
				}
				if (ind != 4*cntSegSeg) Fatal("Group%d arm reduction array size error!", subGroup-GROUP0);
			
			}
			
			// Allocate memory on device
			if (subGroup == GROUP2) {
				device->nSegSeg2 = cntSegSeg;
				device->blocksSegSegs2 = blocksSegSegs;
			
				HANDLE_ERROR(cudaMalloc(&device->f2, sizeof(double3)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g2, sizeof(int2)*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g2pos, sizeof(int2)*nodeCount));
				HANDLE_ERROR(cudaMalloc(&device->g2ind, sizeof(int)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g2arms_pos, sizeof(int2)*armCount));
				HANDLE_ERROR(cudaMalloc(&device->g2arms_ind, sizeof(int)*4*cntSegSeg));
				
				HANDLE_ERROR(cudaMemcpy(device->g2, gList, sizeof(int2)*cntSegSeg, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g2pos, redg_pos, sizeof(int2)*nodeCount, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g2ind, redg_ind, sizeof(int)*4*cntSegSeg, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g2arms_pos, garm_pos, sizeof(int2)*armCount, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g2arms_ind, garm_ind, sizeof(int)*4*cntSegSeg, cudaMemcpyHostToDevice));
			
			} else if (subGroup == GROUP3) {
				device->nSegSeg3 = cntSegSeg;
				device->blocksSegSegs3 = blocksSegSegs;
			
				HANDLE_ERROR(cudaMalloc(&device->f3, sizeof(double3)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g3, sizeof(int2)*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g3pos, sizeof(int2)*nodeCount));
				HANDLE_ERROR(cudaMalloc(&device->g3ind, sizeof(int)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g3arms_pos, sizeof(int2)*armCount));
				HANDLE_ERROR(cudaMalloc(&device->g3arms_ind, sizeof(int)*4*cntSegSeg));
				
				HANDLE_ERROR(cudaMemcpy(device->g3, gList, sizeof(int2)*cntSegSeg, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g3pos, redg_pos, sizeof(int2)*nodeCount, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g3ind, redg_ind, sizeof(int)*4*cntSegSeg, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g3arms_pos, garm_pos, sizeof(int2)*armCount, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g3arms_ind, garm_ind, sizeof(int)*4*cntSegSeg, cudaMemcpyHostToDevice));
			
			} else if (subGroup == GROUP4) {
				device->nSegSeg4 = cntSegSeg;
				device->blocksSegSegs4 = blocksSegSegs;
				
				HANDLE_ERROR(cudaMalloc(&device->f4, sizeof(double3)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g4, sizeof(int2)*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g4pos, sizeof(int2)*nodeCount));
				HANDLE_ERROR(cudaMalloc(&device->g4ind, sizeof(int)*4*cntSegSeg));
				HANDLE_ERROR(cudaMalloc(&device->g4arms_pos, sizeof(int2)*armCount));
				HANDLE_ERROR(cudaMalloc(&device->g4arms_ind, sizeof(int)*4*cntSegSeg));
				
				HANDLE_ERROR(cudaMemcpy(device->g4, gList, sizeof(int2)*cntSegSeg, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g4pos, redg_pos, sizeof(int2)*nodeCount, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g4ind, redg_ind, sizeof(int)*4*cntSegSeg, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g4arms_pos, garm_pos, sizeof(int2)*armCount, cudaMemcpyHostToDevice));
				HANDLE_ERROR(cudaMemcpy(device->g4arms_ind, garm_ind, sizeof(int)*4*cntSegSeg, cudaMemcpyHostToDevice));
			}
			
			free(gList);
			free(redg_pos);
			free(redg_ind);
			free(garm_pos);
			free(garm_ind);
		//}
		
		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    MoveInteractionGroup
 *      Description: Move a given segment pair interaction from group 0
 *                   to another group
 *
 *-----------------------------------------------------------------------*/
void MoveInteractionGroup(Subcyc_t *subcyc, int SegSegID, int subGroup)
{
		
		if (subGroup == 1) {
			Fatal("Should not move interaction to group %d!", subGroup);
			
		} else if (subGroup == 2) {
			
			if (subcyc->SegSegListG2 == NULL || 
				subcyc->SegSegListG2_cnt >= subcyc->SegSegListG2_siz) {
				subcyc->SegSegListG2_siz += 1000;
				subcyc->SegSegListG2 = (SegSeg_t *)realloc(subcyc->SegSegListG2,
										sizeof(SegSeg_t) * subcyc->SegSegListG2_siz);
			}
			
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].seg1  = subcyc->SegSegListG0[SegSegID].seg1;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].seg2  = subcyc->SegSegListG0[SegSegID].seg2;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].flag  = 1;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].setSeg1Forces = 1;
			subcyc->SegSegListG2[subcyc->SegSegListG2_cnt].setSeg2Forces = 1;
			subcyc->SegSegListG2_cnt++;
			
		} else if (subGroup == 3) {
			
			if (subcyc->SegSegListG3 == NULL || 
				subcyc->SegSegListG3_cnt >= subcyc->SegSegListG3_siz) {
				subcyc->SegSegListG3_siz += 1000;
				subcyc->SegSegListG3 = (SegSeg_t *)realloc(subcyc->SegSegListG3,
										sizeof(SegSeg_t) * subcyc->SegSegListG3_siz);
			}
			
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].seg1  = subcyc->SegSegListG0[SegSegID].seg1;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].seg2  = subcyc->SegSegListG0[SegSegID].seg2;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].flag  = 1;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].setSeg1Forces = 1;
			subcyc->SegSegListG3[subcyc->SegSegListG3_cnt].setSeg2Forces = 1;
			subcyc->SegSegListG3_cnt++;
			
		} else if (subGroup == 4) {
			
			if (subcyc->SegSegListG4 == NULL || 
				subcyc->SegSegListG4_cnt >= subcyc->SegSegListG4_siz) {
				subcyc->SegSegListG4_siz += 1000;
				subcyc->SegSegListG4 = (SegSeg_t *)realloc(subcyc->SegSegListG4,
										sizeof(SegSeg_t) * subcyc->SegSegListG4_siz);
			}
			
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].seg1  = subcyc->SegSegListG0[SegSegID].seg1;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].seg2  = subcyc->SegSegListG0[SegSegID].seg2;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].flag  = 1;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].setSeg1Forces = 1;
			subcyc->SegSegListG4[subcyc->SegSegListG4_cnt].setSeg2Forces = 1;
			subcyc->SegSegListG4_cnt++;
			
		} else {
			Fatal("Unknown subGroup %d!", subGroup);
		}
							
		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    SubcycleIntegratorGPU
 *      Description: Perform sublcycling time-integration on the GPU
 *
 *-----------------------------------------------------------------------*/
void SubcycleIntegratorGPU(Home_t *home)
{
		int       i, j, k, l, n1, n2, mobErr;
		double    mobMatrix[3][3];
		Node_t    *node1, *node2;
		Param_t   *param;
		Subcyc_t  *subcyc;
		
		param = home->param;
		subcyc  = home->subcyc;
		
		
		Device_t *device;
		device = (Device_t*)malloc(sizeof(Device_t));
		
		cudaEvent_t start, stop;
		HANDLE_ERROR(cudaEventCreate(&start));
		HANDLE_ERROR(cudaEventCreate(&stop));
		HANDLE_ERROR(cudaEventRecord(start, 0));
		
/*
 *		Set external stress to GPU
 */		
		HANDLE_ERROR(cudaMalloc(&device->esig, sizeof(double)*6));
		HANDLE_ERROR(cudaMemcpy(device->esig, param->appliedStress, sizeof(double)*6, cudaMemcpyHostToDevice));
		
/*
 *		Pack nodes and segments for GPU
 */		
		int nodeCount = 0;
		int segCount = 0;
		int armCount = 0;
		
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			armCount += node1->numNbrs;
			for (j = 0; j < node1->numNbrs; j++) {
				node2 = GetNeighborNode(home, node1, j);
				if (node2 == (Node_t *)NULL) continue;
				if (OrderNodes(node2, node1) != 1) continue;
				segCount++;
			}
			/* id in the GPU node vector */
			node1->subindex = nodeCount;
			nodeCount++;
		}
		if (armCount != 2*segCount) Fatal("armCount != 2*segCount");
		
		
/*
 * 		Determine the number of blocks and threads on which the
 * 		GPU kernels will be executed. Make sure the number of blocks 
 * 		does not exceed the hard limit on GPU devices with compute 
 * 		capability <= 3.x (maxBlocks = 65535).
 */
		int threadsPerBlock, threadsSize[8], maxThreadsSize, blockSize;
		
		threadsSize[0] = nodeCount;
		threadsSize[1] = segCount;
		threadsSize[2] = armCount;
		threadsSize[3] = subcyc->SegListG1_cnt;
		threadsSize[4] = subcyc->SegSegListG0_cnt;
		threadsSize[5] = subcyc->SegSegListG2_cnt;
		threadsSize[6] = subcyc->SegSegListG3_cnt;
		threadsSize[7] = subcyc->SegSegListG4_cnt;
		maxThreadsSize = *std::max_element(threadsSize, threadsSize+8);
		
		GetThreadsPerBlock(home->deviceProp, maxThreadsSize, &threadsPerBlock, &blockSize);
		
		if (blockSize > home->deviceProp->maxBlocks) {
			Fatal("Max number of CUDA blocks exceeded!");
		} else {
			//printf("GPU: threadsPerBlock = %d, blockSize = %d\n", threadsPerBlock, blockSize);
		}
		
		device->threadsPerBlock = threadsPerBlock;
		int blocksNodes = (nodeCount + threadsPerBlock - 1) / threadsPerBlock;
		int blocksSegs = (segCount + threadsPerBlock - 1) / threadsPerBlock;
		int blocksArms = (armCount + threadsPerBlock - 1) / threadsPerBlock;
		
		
		device->nodeCount = nodeCount;
		device->segCount = segCount;
		device->armCount = armCount;
		device->blocksNodes = blocksNodes;
		device->blocksSegs = blocksSegs;
		device->blocksArms = blocksArms;
		
		double3 *b = (double3*)malloc(segCount*sizeof(double3));
		double3 *r = (double3*)malloc(nodeCount*sizeof(double3));
		double  *mob = (double*)malloc(9*nodeCount*sizeof(double));
		int     *armid = (int*)malloc(sizeof(int)*4*segCount);
		
		double3 *cc = (double3*)malloc(nodeCount*sizeof(double3));
		int     cellX, cellY, cellZ;
        double  xCenter, yCenter, zCenter;
        Cell_t  *cell;
        
		
		nodeCount = 0;
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			r[nodeCount].x = node1->x;
			r[nodeCount].y = node1->y;
			r[nodeCount].z = node1->z;
			
			/* Mobility matrix */
			mobErr = Mobility_FCC_0_matrix_GPU(home, node1, mobMatrix);
			for (k = 0; k < 3; k++)
				for (l = 0; l < 3; l++)
					mob[nodeCount*9+k*3+l] = mobMatrix[k][l];
			
			/* Find cell center */ // WARNING
			cell = home->cellKeys[node1->cellIdx];
			cellX = cell->xIndex;
			cellY = cell->yIndex;
			cellZ = cell->zIndex;
			FindCellCenter(param, (real8)(cellX-1), (real8)(cellY-1),
                          (real8)(cellZ-1), 2, &xCenter, &yCenter, &zCenter);
			cc[nodeCount].x = xCenter;
			cc[nodeCount].y = yCenter;
			cc[nodeCount].z = zCenter;
			/* Find cell center */
			
			nodeCount++;
		}
		
		HANDLE_ERROR(cudaMalloc(&device->cc, sizeof(double3)*nodeCount));
		HANDLE_ERROR(cudaMemcpy(device->cc, cc, sizeof(double3)*nodeCount, cudaMemcpyHostToDevice));
		
/*
 *		Pack remote force information (FMM) for GPU
 */	
		if (param->elasticinteraction && param->forceCutOff == 0 && param->fmEnabled) {
			
			subcyc->sigbFMM = (real8**)malloc(segCount*sizeof(real8*));
			RemoteSigbSub(home);
			
			int numPoints = param->fmNumPoints;
			if (numPoints != 3) Fatal("fmNumPoints needs to be 3 for GPU subcycling");
			
			double3 *fmm = (double3*)malloc(segCount*numPoints*sizeof(double3));
			
			for (i = 0; i < segCount; i++) {
				for (j = 0; j < numPoints; j++) {
					fmm[i*numPoints+j].x = subcyc->sigbFMM[i][j*3+0];
					fmm[i*numPoints+j].y = subcyc->sigbFMM[i][j*3+1];
					fmm[i*numPoints+j].z = subcyc->sigbFMM[i][j*3+2];
				}
			}
			
			HANDLE_ERROR(cudaMalloc(&device->fmm, sizeof(double3)*segCount*numPoints));
			HANDLE_ERROR(cudaMemcpy(device->fmm, fmm, sizeof(double3)*segCount*numPoints, cudaMemcpyHostToDevice));
			
			free(fmm);
			for (i = 0; i < subcyc->SegListG1_cnt; i++) {
				free(subcyc->sigbFMM[i]);
			}
			free(subcyc->sigbFMM);
		}
		
/*
 *		Pack group1 and segments for GPU
 */		
		if (segCount != subcyc->SegListG1_cnt) Fatal("segCount != subcyc->SegListG1_cnt");
		Segm_t *SegList = subcyc->SegListG1;
		
		std::vector<std::vector<int> > redg1(nodeCount);
		int2 *g1List = (int2*)malloc(segCount*sizeof(int2));
		int  *g1arms = (int*)malloc(armCount*sizeof(int));
		
		for (j = 0; j < subcyc->SegListG1_cnt; j++) {
			
			node1 = SegList[j].seg->node1;
			node2 = SegList[j].seg->node2;
			
			n1 = node1->subindex;
			n2 = node2->subindex;
			
			redg1[n1].push_back(j*2);
			redg1[n2].push_back(j*2+1);
			
			g1List[j].x = n1;
			g1List[j].y = n2;
			
			k = SegList[j].seg->armID12;
			l = SegList[j].seg->armID21;
			
			b[j].x = node1->burgX[k];
			b[j].y = node1->burgY[k];
			b[j].z = node1->burgZ[k];
			
			node1->armid[k] = 2*j+0;
			node2->armid[l] = 2*j+1;
			
			g1arms[node1->armid[k]] = j*2;
			g1arms[node2->armid[l]] = j*2+1;
			
			armid[(node1->armid[k])*2+0] = node1->myTag.index;
			armid[(node1->armid[k])*2+1] = k;
			armid[(node2->armid[l])*2+0] = node2->myTag.index;
			armid[(node2->armid[l])*2+1] = l;
			
			/* id in the GPU seg vector */
			SegList[j].seg->subindex = j;
		}
		
		// Pack reduction array
		int2 *redg1_pos = (int2*)malloc(nodeCount*sizeof(int2));
		int  *redg1_ind = (int*)malloc(2*segCount*sizeof(int));
		int ind = 0;
		for (j = 0; j < nodeCount; j++) {
			redg1_pos[j].x = ind;
			for (k = 0; k < redg1[j].size(); k++) {
				redg1_ind[ind++] = redg1[j][k];
			}
			redg1_pos[j].y = ind;
		}
		if (ind != 2*segCount) Fatal("Group1 reduction array size error!");
		
		HANDLE_ERROR(cudaMalloc(&device->g1, sizeof(int2)*segCount));
		HANDLE_ERROR(cudaMalloc(&device->g1pos, sizeof(int2)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->g1ind, sizeof(int)*2*segCount));
		HANDLE_ERROR(cudaMalloc(&device->g1arms, sizeof(int)*armCount));
					
		HANDLE_ERROR(cudaMemcpy(device->g1, g1List, sizeof(int2)*segCount, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(device->g1pos, redg1_pos, sizeof(int2)*nodeCount, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(device->g1ind, redg1_ind, sizeof(int)*2*segCount, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(device->g1arms, g1arms, sizeof(int)*armCount, cudaMemcpyHostToDevice));
		
		free(g1List);
		free(g1arms);
		free(redg1_pos);
		free(redg1_ind);
		
/*
 *		Memory allocation on GPU
 */		
		HANDLE_ERROR(cudaMalloc(&device->r, sizeof(double3)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->r0, sizeof(double3)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->n, sizeof(int)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->b, sizeof(double3)*segCount));
		HANDLE_ERROR(cudaMalloc(&device->mob, sizeof(double)*9*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->fseg, sizeof(double3)*2*segCount));
		HANDLE_ERROR(cudaMalloc(&device->Bseg, sizeof(double)*segCount));
		HANDLE_ERROR(cudaMalloc(&device->f, sizeof(double3)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->B, sizeof(double)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->v, sizeof(double3)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->v0, sizeof(double3)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->rkf, sizeof(double3)*5*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->farms, sizeof(double3)*armCount));
		
		// Error arrays
		HANDLE_ERROR(cudaMalloc(&device->e1, sizeof(double)*blocksNodes));
		HANDLE_ERROR(cudaMalloc(&device->e2, sizeof(double)*blocksNodes));
		device->e1_ptr = thrust::device_pointer_cast(device->e1);
		device->e2_ptr = thrust::device_pointer_cast(device->e2);
		
/*
 *		Copy memory from host to device
 */
		HANDLE_ERROR(cudaMemcpy(device->r, r, sizeof(double3)*nodeCount, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(device->b, b, sizeof(double3)*segCount, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(device->mob, mob, sizeof(double)*9*nodeCount, cudaMemcpyHostToDevice));
		

/*
 *		Perform RKF integration on GPU. First, pack the
 *      segment pairs and send them to the GPU. Then, calculate
 *      the interaction distances to assign each interaction
 *      into a group. After this, time integrate group 0 forces.
 */		
		SendSubgroupGPU(home, device, GROUP0);
		SegSegListGPU(home, device);
		RKFIntegratorGPU(home, device, GROUP0, 0);

/*
 *		Update the different groups by moving segment pairs
 *      interactions that have been flagged.
 */		
		int *g0flag = (int*)malloc(device->nSegSeg0*sizeof(int));
		HANDLE_ERROR(cudaMemcpy(g0flag, device->g0flag, sizeof(int)*device->nSegSeg0, cudaMemcpyDeviceToHost));
		
		int cntSegSeg = 0;
		for (j = 0; j < subcyc->SegSegListG0_cnt; j++) {
			if (subcyc->SegSegListG0[j].flag == 0) continue;
			
			int subGroup = g0flag[cntSegSeg];
			cntSegSeg++;
			
			if (subGroup < 0) {
				/* Interactions that are flagged -1 are to be 
				 * ignored (e.g. very small segments) */
				subcyc->SegSegListG0[j].flag = 0;
				continue;
				
			} else if (subGroup > 0) {
				MoveInteractionGroup(subcyc, j, subGroup);
				subcyc->SegSegListG0[j].flag = 0;
			}
		}
		free(g0flag);
			
		SendSubgroupGPU(home, device, GROUP4);
		SendSubgroupGPU(home, device, GROUP3);
		SendSubgroupGPU(home, device, GROUP2);
		
/*
 *		Time integrate group 1, 2, 3 and 4 interactions (subcycle).
 */		
		real8   subTime1, subTime2, subTime3, subTime4;
        real8   totalsubDT, nextDTsub, oldDTsub, newDTsub;
        int     subGroup, cutDT;
        
        //Initialize the time for each group based on whether it has any forces in it
        if (segCount > 0) subTime1 = 0.0;
        else              subTime1 = param->realdt;
		
        if (device->nSegSeg2 > 0) subTime2 = 0.0;
        else                      subTime2 = param->realdt;
		
        if (device->nSegSeg3 > 0) subTime3 = 0.0;
        else                      subTime3 = param->realdt;
		
        if (device->nSegSeg4 > 0) subTime4 = 0.0;
        else                      subTime4 = param->realdt;
		
        //Initialize some other stuff
		if (home->cycle == 0) nextDTsub = param->realdt;
		totalsubDT = 0.0;
        subcyc->numSubCycle1 = 0;
        subcyc->numSubCycle2 = 0;
        subcyc->numSubCycle3 = 0;
        subcyc->numSubCycle4 = 0;
        int oldGroup = -1;
        int nSubcyc;
        int totSubcyc = 0;

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

            //If we switched groups, reset subcycle count
            if (subGroup != oldGroup) nSubcyc = 0;
            oldGroup = subGroup;

            //Make sure we don't pass the global group in time
			if (totalsubDT + nextDTsub > param->realdt) {
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
			RKFIntegratorGPU(home, device, subGroup, nSubcyc);
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
			
			totSubcyc++;
		}
		
	
/*
 *		We are done with subcycling. Now recalculate all the
 *      forces and mobilities.
 */		
		RKFIntegratorGPU(home, device, FULL, -1);

/*
 *		Unpack arms forces back on the CPU
 */		
		double3 *farms = (double3*)malloc(sizeof(double3)*armCount);
		HANDLE_ERROR(cudaMemcpy(farms, device->farms, sizeof(double3)*armCount, cudaMemcpyDeviceToHost));
		
		for (i = 0; i < armCount; i++) {
			j = armid[2*i+0];
			k = armid[2*i+1];
			node1 = home->nodeKeys[j];
			node1->armfx[k] = farms[i].x;
			node1->armfy[k] = farms[i].y;
			node1->armfz[k] = farms[i].z;
		}
		
		free(farms);
		
/*
 *		Copy new nodal positions, forces, and velocities back to the CPU
 */	
		HANDLE_ERROR(cudaMemcpy(r, device->r, sizeof(double3)*nodeCount, cudaMemcpyDeviceToHost));
		
		double3 *f = (double3*)malloc(sizeof(double3)*nodeCount);
		double3 *v = (double3*)malloc(sizeof(double3)*nodeCount);
		HANDLE_ERROR(cudaMemcpy(f, device->f, sizeof(double3)*nodeCount, cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(v, device->v, sizeof(double3)*nodeCount, cudaMemcpyDeviceToHost));
		
		nodeCount = 0;
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			node1->x = r[nodeCount].x;
			node1->y = r[nodeCount].y;
			node1->z = r[nodeCount].z;
			
			node1->fX = f[nodeCount].x;
			node1->fY = f[nodeCount].y;
			node1->fZ = f[nodeCount].z;
			
			node1->vX = v[nodeCount].x;
			node1->vY = v[nodeCount].y;
			node1->vZ = v[nodeCount].z;
			
			node1->flags &= (~NODE_RESET_FORCES);
			
			nodeCount++;
		}
		
		
/*
 *		Check that no errors have been encountered. Errors can occur
 *      when the GPU memory becomes full or fragmented. When this
 *      happens we must stop the simulation, otherwise all subsequent
 *      results will be garbage.
 */		
		CheckErrorGPU("SubcycleIntegratorGPU");
		
/*
 *		Free memory
 */		
		HANDLE_ERROR(cudaFree(device->r));
		HANDLE_ERROR(cudaFree(device->r0));
		HANDLE_ERROR(cudaFree(device->n));
		HANDLE_ERROR(cudaFree(device->b));
		HANDLE_ERROR(cudaFree(device->mob));
		HANDLE_ERROR(cudaFree(device->fseg));
		HANDLE_ERROR(cudaFree(device->Bseg));
		HANDLE_ERROR(cudaFree(device->f));
		HANDLE_ERROR(cudaFree(device->B));
		HANDLE_ERROR(cudaFree(device->v));
		HANDLE_ERROR(cudaFree(device->v0));
		HANDLE_ERROR(cudaFree(device->rkf));
		
		HANDLE_ERROR(cudaFree(device->esig));
		
		HANDLE_ERROR(cudaFree(device->g1));
		HANDLE_ERROR(cudaFree(device->g1pos));
		HANDLE_ERROR(cudaFree(device->g1ind));
		HANDLE_ERROR(cudaFree(device->g1arms));
		
		HANDLE_ERROR(cudaFree(device->farms));
		
		//if (param->elasticinteraction) {
			HANDLE_ERROR(cudaFree(device->f0));
			HANDLE_ERROR(cudaFree(device->g0));
			HANDLE_ERROR(cudaFree(device->g0pos));
			HANDLE_ERROR(cudaFree(device->g0ind));
			HANDLE_ERROR(cudaFree(device->g0flag));
			HANDLE_ERROR(cudaFree(device->g0dist2));
			HANDLE_ERROR(cudaFree(device->g0arms_pos));
			HANDLE_ERROR(cudaFree(device->g0arms_ind));
			
			HANDLE_ERROR(cudaFree(device->f2));
			HANDLE_ERROR(cudaFree(device->g2));
			HANDLE_ERROR(cudaFree(device->g2pos));
			HANDLE_ERROR(cudaFree(device->g2ind));
			HANDLE_ERROR(cudaFree(device->g2arms_pos));
			HANDLE_ERROR(cudaFree(device->g2arms_ind));
			
			HANDLE_ERROR(cudaFree(device->f3));
			HANDLE_ERROR(cudaFree(device->g3));
			HANDLE_ERROR(cudaFree(device->g3pos));
			HANDLE_ERROR(cudaFree(device->g3ind));
			HANDLE_ERROR(cudaFree(device->g3arms_pos));
			HANDLE_ERROR(cudaFree(device->g3arms_ind));
			
			HANDLE_ERROR(cudaFree(device->f4));
			HANDLE_ERROR(cudaFree(device->g4));
			HANDLE_ERROR(cudaFree(device->g4pos));
			HANDLE_ERROR(cudaFree(device->g4ind));
			HANDLE_ERROR(cudaFree(device->g4arms_pos));
			HANDLE_ERROR(cudaFree(device->g4arms_ind));
		//}
			
		if (param->elasticinteraction) {	
			if (param->forceCutOff == 0 && param->fmEnabled) {
				HANDLE_ERROR(cudaFree(device->fmm));
			}
		}
		
		HANDLE_ERROR(cudaFree(device->e1));
		HANDLE_ERROR(cudaFree(device->e2));
		
		free(r);
		free(b);
		free(mob);
		free(f);
		free(v);
		
		HANDLE_ERROR(cudaFree(device->cc));
		free(cc);
		
		free(armid);
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			free(node1->armid);
			free(node1->armInt);
		}
		
		free(device);
		
		HANDLE_ERROR(cudaEventRecord(stop, 0));
		HANDLE_ERROR(cudaEventSynchronize(stop));
		float gputime;
		HANDLE_ERROR(cudaEventElapsedTime(&gputime, start, stop));
		//printf("SubcycleIntegratorGPU time: %f ms\n", gputime);
		
		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    PairSegSegForcesGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void PairSegSegForcesGPU(int pairCount, double3 *r, int *n12, int *n34, double *cc, double *b1, double *b2, double *fpair)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < pairCount) {
		
		int arm = n34[i*3+0];
		int node3 = n34[i*3+1];
		int node4 = n34[i*3+2];
		
		int node1 = n12[arm*2+0];
		int node2 = n12[arm*2+1];
		
		double3 r1, r2, r3, r4;
		double3 bs1, bs2;
		double3 f1, f2, f3, f4;
		
		bs1.x = b1[arm*3+0];
		bs1.y = b1[arm*3+1];
		bs1.z = b1[arm*3+2];
		
		bs2.x = b2[i*3+0];
		bs2.y = b2[i*3+1];
		bs2.z = b2[i*3+2];
		
		r1 = r[node1];
		r2 = r[node2];
		r3 = r[node3];
		r4 = r[node4];
		
		// Cell center
		double3 rc;
		rc.x = cc[arm*3+0];
		rc.y = cc[arm*3+1];
		rc.z = cc[arm*3+2];
		
		double3 dr1, dr2;
		double  L1s, L2s;
		dr1.x = r2.x - r1.x;
		dr1.y = r2.y - r1.y;
		dr1.z = r2.z - r1.z;
		ZImageGPU(dr1);
		L1s = dr1.x*dr1.x + dr1.y*dr1.y + dr1.z*dr1.z;
		
		dr2.x = r4.x - r3.x;
		dr2.y = r4.y - r3.y;
		dr2.z = r4.z - r3.z;
		ZImageGPU(dr2);
		L2s = dr2.x*dr2.x + dr2.y*dr2.y + dr2.z*dr2.z;
		
		if (L1s < 1.e-20 || L2s < 1.e-20) {
			
			f1.x = 0.0; f1.y = 0.0; f1.z = 0.0;
			f2.x = 0.0; f2.y = 0.0; f2.z = 0.0;
			
		} else {
		
			r2.x = r1.x + dr1.x;
			r2.y = r1.y + dr1.y;
			r2.z = r1.z + dr1.z;
			
			// Cell center
			dr1.x = r3.x - rc.x;
			dr1.y = r3.y - rc.y;
			dr1.z = r3.z - rc.z;
			ZImageGPU(dr1);
			
			r3.x = rc.x + dr1.x;
			r3.y = rc.y + dr1.y;
			r3.z = rc.z + dr1.z;
			
			dr2.x = r4.x - r3.x;
			dr2.y = r4.y - r3.y;
			dr2.z = r4.z - r3.z;
			ZImageGPU(dr2);
			
			r4.x = r3.x + dr2.x;
			r4.y = r3.y + dr2.y;
			r4.z = r3.z + dr2.z;
			
			SegSegForceIsotropicGPU(r1, r2, r3, r4, bs1, bs2, f1, f2, f3, f4);
			
		}
		
		fpair[i*6+0] = f1.x;
		fpair[i*6+1] = f1.y;
		fpair[i*6+2] = f1.z;
		fpair[i*6+3] = f2.x;
		fpair[i*6+4] = f2.y;
		fpair[i*6+5] = f2.z;
		
	}
	
	return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    DevicePair_t
 *
 *-----------------------------------------------------------------------*/
typedef struct {
	int      nodeSize, pairSize, armSize;
	double3  *r_host, *r;
	double   *fpair;
	int      *n12, *n34;
	double   *cc, *b1, *b2;
} DevicePair_t;

DevicePair_t *devicePair;

/*------------------------------------------------------------------------
 *
 *      Function:    InitializeNodeForceGPU
 *
 *-----------------------------------------------------------------------*/
void InitializeNodeForceGPU(Home_t *home)
{
	devicePair = (DevicePair_t*)malloc(sizeof(DevicePair_t));
	devicePair->nodeSize = 0;
	devicePair->pairSize = 0;
	devicePair->armSize = 0;
	
	return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    FinalizeNodeForceGPU
 *
 *-----------------------------------------------------------------------*/
void FinalizeNodeForceGPU(Home_t *home)
{
	if (devicePair->nodeSize > 0) {
		free(devicePair->r_host);
		HANDLE_ERROR(cudaFree(devicePair->r));
	}
	if (devicePair->pairSize > 0) {
		HANDLE_ERROR(cudaFree(devicePair->fpair));
		HANDLE_ERROR(cudaFree(devicePair->n34));
		HANDLE_ERROR(cudaFree(devicePair->b2));
	}
	if (devicePair->armSize > 0) {
		HANDLE_ERROR(cudaFree(devicePair->cc));
		HANDLE_ERROR(cudaFree(devicePair->b1));
		HANDLE_ERROR(cudaFree(devicePair->n12));
	}
	free(devicePair);
	
	return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    SetOneNodeForceGPU
 *      Description: Compute nodal forces on nodes being splitted.
 *                   Note: the implementation has been optimized to 
 *                   reduce memory overheads.
 *
 *-----------------------------------------------------------------------*/
void SetOneNodeForceGPU(Home_t *home, SplitSegSeg_t *splitSegSegList)
{
		int     i, nodeSize, pairSize, armSize;
		int     segPairCnt, nodeCount;
		Node_t  *node1;
		
		segPairCnt = splitSegSegList->segPairCnt;
		nodeCount = splitSegSegList->nodeCount;
		
		// Pack the nodes
		nodeSize = devicePair->nodeSize;
		if (nodeCount > nodeSize) {
			if (nodeSize > 0) {
				free(devicePair->r_host);
				HANDLE_ERROR(cudaFree(devicePair->r));
			}
			nodeSize = 2*nodeCount;
			devicePair->nodeSize = nodeSize;
			devicePair->r_host = (double3*)malloc(sizeof(double3)*nodeSize);
			HANDLE_ERROR(cudaMalloc(&devicePair->r, sizeof(double3)*nodeSize));
		}
		
		nodeCount = 0;
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			devicePair->r_host[nodeCount].x = node1->x;
			devicePair->r_host[nodeCount].y = node1->y;
			devicePair->r_host[nodeCount].z = node1->z;
			nodeCount++;
		}
		
		HANDLE_ERROR(cudaMemcpy(devicePair->r, devicePair->r_host, sizeof(double3)*nodeCount, cudaMemcpyHostToDevice));
		
		
		// Pack the interactions
		pairSize = devicePair->pairSize;
		if (segPairCnt > pairSize) {
			if (pairSize > 0) {
				HANDLE_ERROR(cudaFree(devicePair->fpair));
				HANDLE_ERROR(cudaFree(devicePair->n34));
				HANDLE_ERROR(cudaFree(devicePair->b2));
			}
			pairSize = 2*segPairCnt;
			devicePair->pairSize = pairSize;
			HANDLE_ERROR(cudaMalloc(&devicePair->fpair, sizeof(double)*6*pairSize));
			HANDLE_ERROR(cudaMalloc(&devicePair->n34, sizeof(int)*pairSize*3));
			HANDLE_ERROR(cudaMalloc(&devicePair->b2, sizeof(double)*pairSize*3));
		}
		
		int nArms = splitSegSegList->nArms;
		armSize = devicePair->armSize;
		if (nArms > armSize) {
			if (armSize > 0) {
				HANDLE_ERROR(cudaFree(devicePair->cc));
				HANDLE_ERROR(cudaFree(devicePair->b1));
				HANDLE_ERROR(cudaFree(devicePair->n12));
			}
			armSize = 2*nArms;
			devicePair->armSize = armSize;
			HANDLE_ERROR(cudaMalloc(&devicePair->cc, sizeof(double)*armSize*3));
			HANDLE_ERROR(cudaMalloc(&devicePair->b1, sizeof(double)*armSize*3));
			HANDLE_ERROR(cudaMalloc(&devicePair->n12, sizeof(int)*armSize*2));
		}
		
		HANDLE_ERROR(cudaMemcpy(devicePair->cc, splitSegSegList->cc, sizeof(double)*nArms*3, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(devicePair->b1, splitSegSegList->b1, sizeof(double)*nArms*3, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(devicePair->n12, splitSegSegList->n12, sizeof(int)*nArms*2, cudaMemcpyHostToDevice));
		
		HANDLE_ERROR(cudaMemcpy(devicePair->n34, splitSegSegList->n34, sizeof(int)*segPairCnt*3, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(devicePair->b2, splitSegSegList->b2, sizeof(double)*segPairCnt*3, cudaMemcpyHostToDevice));
			
		
		// Compute segment / segment interactions
		int threadsPerBlock, blocksPairs;
		GetThreadsPerBlock(home->deviceProp, segPairCnt, &threadsPerBlock, &blocksPairs);
		
		PairSegSegForcesGPU<<<blocksPairs,threadsPerBlock>>>(segPairCnt, devicePair->r, devicePair->n12, devicePair->n34, devicePair->cc, devicePair->b1, devicePair->b2, devicePair->fpair);
		
		HANDLE_ERROR(cudaMemcpy(splitSegSegList->fpair, devicePair->fpair, sizeof(double)*6*segPairCnt, cudaMemcpyDeviceToHost));
		
		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    LocalSegSegForcesGPU
 *
 *-----------------------------------------------------------------------*/
__global__ void LocalSegSegForcesGPU(int pairCount, double3 *r, int4 *pair, double3 *b1, double3 *b2, double3 *cc, double *fpair)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < pairCount) {
		
		int4 nodes = pair[i];
		
		double3 r1, r2, r3, r4;
		double3 bs1, bs2;
		double3 f1, f2, f3, f4;
		
		bs1 = b1[i];
		bs2 = b2[i];
		
		r1 = r[nodes.x];
		r2 = r[nodes.y];
		r3 = r[nodes.z];
		r4 = r[nodes.w];
		
		double3 dr1, dr2;
		double  L1s, L2s;
		dr1.x = r2.x - r1.x;
		dr1.y = r2.y - r1.y;
		dr1.z = r2.z - r1.z;
		ZImageGPU(dr1);
		L1s = dr1.x*dr1.x + dr1.y*dr1.y + dr1.z*dr1.z;
		
		dr2.x = r4.x - r3.x;
		dr2.y = r4.y - r3.y;
		dr2.z = r4.z - r3.z;
		ZImageGPU(dr2);
		L2s = dr2.x*dr2.x + dr2.y*dr2.y + dr2.z*dr2.z;
		
		if (L1s < 1.e-20 || L2s < 1.e-20) {
			
			f1.x = 0.0; f1.y = 0.0; f1.z = 0.0;
			f2.x = 0.0; f2.y = 0.0; f2.z = 0.0;
			f3.x = 0.0; f3.y = 0.0; f3.z = 0.0;
			f4.x = 0.0; f4.y = 0.0; f4.z = 0.0;
			
		} else {
		
			r2.x = r1.x + dr1.x;
			r2.y = r1.y + dr1.y;
			r2.z = r1.z + dr1.z;
			
			// Cell center
			double3 rc;
			rc = cc[nodes.x];
			dr1.x = r3.x - rc.x;
			dr1.y = r3.y - rc.y;
			dr1.z = r3.z - rc.z;
			ZImageGPU(dr1);
			
			r3.x = rc.x + dr1.x;
			r3.y = rc.y + dr1.y;
			r3.z = rc.z + dr1.z;
			
			dr2.x = r4.x - r3.x;
			dr2.y = r4.y - r3.y;
			dr2.z = r4.z - r3.z;
			ZImageGPU(dr2);
			
			r4.x = r3.x + dr2.x;
			r4.y = r3.y + dr2.y;
			r4.z = r3.z + dr2.z;
			
			SegSegForceIsotropicGPU(r1, r2, r3, r4, bs1, bs2, f1, f2, f3, f4);
			
		}
		
		fpair[i*12+0] = f1.x;
		fpair[i*12+1] = f1.y;
		fpair[i*12+2] = f1.z;
		fpair[i*12+3] = f2.x;
		fpair[i*12+4] = f2.y;
		fpair[i*12+5] = f2.z;
		fpair[i*12+6] = f3.x;
		fpair[i*12+7] = f3.y;
		fpair[i*12+8] = f3.z;
		fpair[i*12+9] = f4.x;
		fpair[i*12+10] = f4.y;
		fpair[i*12+11] = f4.z;
		
	}
	
	return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    LocalSegForcesGPU
 *      Description: Compute the segment pairs forces required in the
 *                   LocalSegForces calculation using the GPU device.
 *
 *-----------------------------------------------------------------------*/
void LocalSegForcesGPU(Home_t *home, int segPairListCnt, SegmentPair_t *segPairList, double *fpair)
{
		int     i, armID12, armID34;
		int     cellX, cellY, cellZ;
		double  xCenter, yCenter, zCenter;
		Node_t  *node1, *node2, *node3, *node4;
		Cell_t  *cell;
		Param_t *param;
		
		param = home->param;
		
		if (segPairListCnt == 0) {
			return;
		}
		
		// Pack nodes
		int nodeCount = 0;
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node1->subindex = nodeCount;
			nodeCount++;
		}
		
		double3 *r = (double3*)malloc(sizeof(double3)*nodeCount);
		double3 *cc = (double3*)malloc(sizeof(double3)*nodeCount);
		
		nodeCount = 0;
		for (i = 0; i < home->newNodeKeyPtr; i++) {
			if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			r[nodeCount].x = node1->x;
			r[nodeCount].y = node1->y;
			r[nodeCount].z = node1->z;
			
			/* Find cell center */
			cell = home->cellKeys[node1->cellIdx];
			cellX = cell->xIndex;
			cellY = cell->yIndex;
			cellZ = cell->zIndex;
			FindCellCenter(param, (real8)(cellX-1), (real8)(cellY-1),
                          (real8)(cellZ-1), 2, &xCenter, &yCenter, &zCenter);
			cc[nodeCount].x = xCenter;
			cc[nodeCount].y = yCenter;
			cc[nodeCount].z = zCenter;
			
			nodeCount++;
		}
		
		Device_t *device;
		device = (Device_t*)malloc(sizeof(Device_t));
		
		HANDLE_ERROR(cudaMalloc(&device->r, sizeof(double3)*nodeCount));
		HANDLE_ERROR(cudaMalloc(&device->cc, sizeof(double3)*nodeCount));
		
		HANDLE_ERROR(cudaMemcpy(device->r, r, sizeof(double3)*nodeCount, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(device->cc, cc, sizeof(double3)*nodeCount, cudaMemcpyHostToDevice));
		
		
		// Pack segment / segment interactions
		int4 *pair = (int4*)malloc(sizeof(int4)*segPairListCnt);
		double3 *b1 = (double3*)malloc(sizeof(double3)*segPairListCnt);
		double3 *b2 = (double3*)malloc(sizeof(double3)*segPairListCnt);
		
		for (i = 0; i < segPairListCnt; i++) {
			node1 = segPairList[i].seg1->node1;
			node2 = segPairList[i].seg1->node2;
			node3 = segPairList[i].seg2->node1;
			node4 = segPairList[i].seg2->node2;
			
			pair[i].x = node1->subindex;
			pair[i].y = node2->subindex;
			pair[i].z = node3->subindex;
			pair[i].w = node4->subindex;
			
			armID12 = GetArmID(home, node1, node2);
			armID34 = GetArmID(home, node3, node4);
			
			b1[i].x = node1->burgX[armID12];
			b1[i].y = node1->burgY[armID12];
			b1[i].z = node1->burgZ[armID12];
			
			b2[i].x = node3->burgX[armID34];
			b2[i].y = node3->burgY[armID34];
			b2[i].z = node3->burgZ[armID34];
		}
		
		HANDLE_ERROR(cudaMalloc(&device->pair, sizeof(int4)*segPairListCnt));
		HANDLE_ERROR(cudaMalloc(&device->b1, sizeof(double3)*segPairListCnt));
		HANDLE_ERROR(cudaMalloc(&device->b2, sizeof(double3)*segPairListCnt));
		
		HANDLE_ERROR(cudaMemcpy(device->pair, pair, sizeof(int4)*segPairListCnt, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(device->b1, b1, sizeof(double3)*segPairListCnt, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(device->b2, b2, sizeof(double3)*segPairListCnt, cudaMemcpyHostToDevice));
		
		
		// Compute segment / segment interactions
		int threadsPerBlock, blocksSegSegs;
		GetThreadsPerBlock(home->deviceProp, segPairListCnt, &threadsPerBlock, &blocksSegSegs);
		
		HANDLE_ERROR(cudaMalloc(&device->fpair, sizeof(double)*12*segPairListCnt));
		
		LocalSegSegForcesGPU<<<blocksSegSegs,threadsPerBlock>>>(segPairListCnt, device->r, device->pair, device->b1, device->b2, device->cc, device->fpair);
		CheckErrorGPU("LocalSegSegForcesGPU");
		
		// Copy back nodal forces
		HANDLE_ERROR(cudaMemcpy(fpair, device->fpair, sizeof(double)*12*segPairListCnt, cudaMemcpyDeviceToHost));
		
		// Free memory
		HANDLE_ERROR(cudaFree(device->r));
		HANDLE_ERROR(cudaFree(device->cc));
		HANDLE_ERROR(cudaFree(device->pair));
		HANDLE_ERROR(cudaFree(device->b1));
		HANDLE_ERROR(cudaFree(device->b2));
		HANDLE_ERROR(cudaFree(device->fpair));
		
		free(r);
		free(cc);
		free(b1);
		free(b2);
		free(pair);
		
		free(device);
		
		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    InitializeParadisGPU
 *      Description: Select the GPU device and set some parameters on it
 *
 *-----------------------------------------------------------------------*/
void InitializeParadisGPU(Home_t *home)
{
	double Lx, Ly, Lz, xc, yc, zc;
	Param_t *param;
	
	param = home->param;
	
	/* Leave if we are not using the GPU subcycling integrator */
	if (strcmp(param->timestepIntegrator, "forceBsubcycle") != 0 || 
	    strcmp(param->subInteg0Integ1, "GPU") != 0) {
		return;
	}
	
#if PARALLEL
	Fatal("GPU subcycling can only be used in serial mode!");
#endif
	if (param->useLabFrame) {
		Fatal("GPU subcycling cannot be used with useLabFrame!");
	}
	if (param->elasticinteraction && param->fmEnabled == 0) {
		Fatal("GPU subcycling cannot be used with Rijm table!");
	}
	if (param->rg1 > 0) {
		Fatal("GPU subcycling cannot be used with rg1 > 0!");
	}
	if (param->mobilityType != MOB_FCC_0) {
		Fatal("GPU subcycling can only be used with MobilityLaw_FCC_0!");
	}
	
	printf("Initializing ParaDiS GPU\n");
	
	DeviceProp_t *deviceProp;
	deviceProp = (DeviceProp_t*)malloc(sizeof(DeviceProp_t));
	SelectCudaDevice(home->deviceID, deviceProp);
	home->deviceProp = deviceProp;
			
	if (param->xBoundType == Periodic) Lx = param->Lx; else Lx = 0.0;
	if (param->yBoundType == Periodic) Ly = param->Ly; else Ly = 0.0;
	if (param->zBoundType == Periodic) Lz = param->Lz; else Lz = 0.0;
		
	xc = (param->maxSideX + param->minSideX) * 0.5;
	yc = (param->maxSideY + param->minSideY) * 0.5;
	zc = (param->maxSideZ + param->minSideZ) * 0.5;
		
	SetVariablesGPU<<<1,1>>>(param->shearModulus, param->pois, param->rc, param->TensionFactor, 
			                 param->MobEdge, param->MobScrew, Lx, Ly, Lz, xc, yc, zc, 
			                 param->rTol, param->rTolth, param->rTolrel, param->Ecore);
	
	return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    GetThreadsPerBlock
 *      Description: Determine the optimal number (lowest) of threads per
 *                   block as a function of the maximum number of blocks
 *
 *-----------------------------------------------------------------------*/
void GetThreadsPerBlock(DeviceProp_t *deviceProp, int threadsSize, int *threadsPerBlock, int *blockSize)
{
	*threadsPerBlock = 64;
	
	while (*threadsPerBlock <= deviceProp->maxThreadsPerBlock) {
		*blockSize = (threadsSize + *threadsPerBlock - 1) / *threadsPerBlock;
		if (*blockSize <= deviceProp->maxBlocks) {
			break;
		}
		*threadsPerBlock *= 2;
	}
		
	return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    GetCudaCores
 *
 *-----------------------------------------------------------------------*/
int GetCudaCores(cudaDeviceProp devProp)
{  
	int cores = 0;
	int mp = devProp.multiProcessorCount;
	switch (devProp.major){
		case 2: // Fermi
			if (devProp.minor == 1) cores = mp * 48;
			else cores = mp * 32;
			break;
		case 3: // Kepler
			cores = mp * 192;
			break;
		case 5: // Maxwell
			cores = mp * 128;
			break;
		case 6: // Pascal
			if (devProp.minor == 1) cores = mp * 128;
			else if (devProp.minor == 0) cores = mp * 64;
			//else printf("Unknown device type\n");
			break;
		default:
			//printf("Unknown device type\n"); 
			break;
	}
	return cores;
}

/*------------------------------------------------------------------------
 *
 *      Function:    SelectCudaDevice
 * 		Description: Select the GPU device
 *
 *-----------------------------------------------------------------------*/
void SelectCudaDevice(int deviceID, DeviceProp_t *deviceProp)
{
	int nDevices, device;
	cudaGetDeviceCount(&nDevices);
	
	printf("\nAvailable GPU device(s): %d\n", nDevices);
	if (nDevices == 0) {
		Fatal("No GPU device is available on this system");
	} else {
		for (int i = 0; i < nDevices; i++) {
			cudaDeviceProp prop;
			cudaGetDeviceProperties(&prop, i);
			printf("  GPU Device ID %d: %s\n", i, prop.name);
		}
	}
	
	if (deviceID == -1) {
		if (nDevices > 1) {
/*
 * 			Loop over the GPU devices and select
 * 			that with the maximum number of processors
 */
			int max_mp = 0;
			for (int i = 0; i < nDevices; i++) {
				cudaDeviceProp prop;
				cudaGetDeviceProperties(&prop, i);
				if (max_mp < prop.multiProcessorCount) {
					max_mp = prop.multiProcessorCount;
					device = i;
				}
			}
			cudaSetDevice(device);
		} else {
			device = 0;
		}
	} else {
/*
 * 		Select the GPU device requested with the -g option
 */
		device = deviceID;
		if (device < 0 || device >= nDevices) {
			printf("\nError: GPU device ID %d is not available on this system\n", device);
			Fatal("Please select a valid device ID in the above list or do not use -g option");
		}
		cudaSetDevice(device);
	}
	
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, device);
	printf("\n**************************************************\n");
	printf("GPU Device ID: %d\n", device);
	printf("  Device name: %s\n", prop.name);
	printf("  Device PCI Bus id: %d\n", prop.pciBusID);
	int cores = GetCudaCores(prop);
	if (cores == 0) {
		printf("  Number of cores: unknown\n");
	} else {
		printf("  Number of cores: %d\n", cores);
	}
	printf("  Clock rate (MHz): %f\n", 1.0*prop.clockRate/1000);
	printf("  Global memory (MB): %f\n", 1.0*prop.totalGlobalMem/1.0e6);
	printf("  Peak Memory Bandwidth (GB/s): %f\n",2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
	printf("  Maximum grid size (x,y,z): %d %d %d\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
	printf("  Maximum threads per block: %d\n", prop.maxThreadsPerBlock);
	printf("**************************************************\n\n");
	
	deviceProp->deviceID = device;
	deviceProp->maxBlocks = prop.maxGridSize[0];
	deviceProp->maxThreadsPerBlock = prop.maxThreadsPerBlock;
}

#endif
