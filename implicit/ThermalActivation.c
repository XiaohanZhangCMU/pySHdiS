/****************************************************************************
 *
 *      Module:         ThermalActivation.c
 *
 *      Author:         Ryan Sills
 *
 *      Description:    This module contains functions for determining
 *						whether a thermally activated event will occur.
 *						This includes calculation of the energy barrier
 *						and event probability, and then using a random
 *						number generator to determine whether the event 
 *						occurs.  
 *
 *
 *      Includes public functions:
 *          ThermalActivation
 *			ParallelTest
 *			FindEscaigDir
 *			CrossSlipStress
 *	    	EnergyBarrierCrossSlipFCC
 *	    
 *
 ***************************************************************************/
#include <math.h>
#include "Home.h"
#include "Util.h"
#include <stdlib.h>

/*---------------------------------------------------------------------------
 *
 *      Function:       ThermalActivation
 *
 *      Description:    Given an energy barrier  and an attempt frequency,
 *						this function evaluates whether a thermally 
 *						activated event occurs.  A 0 or 1 is returned via 
 *						the "doIt" argument. The random number generator has 
 *						already	been initialized in Main.
 *
 *-------------------------------------------------------------------------*/

int ThermalActivation(Param_t *param, real8 Eb, real8 attfreq){

	real8 prob, randnum, dt, T, kB;
	int doIt;	

	dt = param->realdt;
	T = param->TempK;
	kB = 8.6173324e-5;  // eV/K

	prob = dt*attfreq*exp(-Eb/kB/T);

	randnum = (real8)rand() / (real8)RAND_MAX ;

	if (randnum<prob){
		doIt = 1;
	}else{
		doIt = 0;
	}
	return(doIt);

}

/*---------------------------------------------------------------------------
 *
 *      Function:		ParallelTest
 *
 *      Description:    Test if two given vectors are parallel, and return a
 *						logical yes or no (1 or 0). Does not discriminate 
 *						between vectors	pointing in the same or opposite 
 *						directions.
 *
 *-------------------------------------------------------------------------*/

int ParallelTest(real8 vec1[3], real8 vec2[3]) {

	real8 	len1, len2;
	real8	tol = 1e-8;
	int 	test;

    NormalizeVec(vec1);
    NormalizeVec(vec2);

	test = fabs(1.0-fabs(DotProduct(vec1,vec2)))<tol;

	return(test);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    	ThompsonInwardNormal
 *
 *      Description: 	Given a glide plane normal for an FCC metal, return
 *                      the unit normal that points into the Thompson
 *                      tetrahedron.	
 *
 *-------------------------------------------------------------------------*/

void ThompsonInwardPlane(real8 n[3]) {

    real8   sign;
/*
 *  These are the hard-coded glide plane normal directions pointing into the 
 *  Thompson tetrhedron.
 */
    static real8 n1[3] = {1, 1, 1};
    static real8 n2[3] = {1, -1, -1};
    static real8 n3[3] = {-1, 1, -1};
    static real8 n4[3] = {-1, -1, 1};

    if (ParallelTest(n,n1)) {
        sign = (real8) Sign(DotProduct(n,n1));
    } else if (ParallelTest(n,n2)) {
        sign = (real8) Sign(DotProduct(n,n2));
    } else if (ParallelTest(n,n3)) {
        sign = (real8) Sign(DotProduct(n,n3));
    } else if (ParallelTest(n,n4)) {
        sign = (real8) Sign(DotProduct(n,n4));
    } else {
        printf("Warning: Non-(111) plane used for cross-slip " 
		"in ThompsonInwardPlane: %e %e %e\n",
		n[X],n[Y],n[Z]);
        sign = 1.0;
    }    

    n[X] = sign*n[X];
    n[Y] = sign*n[Y];
    n[Z] = sign*n[Z];

    return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    	FindEscaigDir
 *
 *      Description: 	Given a glide plane normal and a Burgers vector, return
 *					 	the positive Escaig stress direction; this is the stress
 *					 	direction that results in widening of the stacking fault.
 *					 	This always corresponds to the direction pointing from the
 *					 	midpoint of the face of the Thomson tetrahedron towards its
 *					 	edge.	
 *
 *-------------------------------------------------------------------------*/

void FindEscaigDir(real8 n[3], real8 b[3], real8 dir[3]) {

	real8	b1[3], b2[3], b3[3];
	real8	dir1[3], dir2[3], dir3[3];

    static real8 n1[3] = {1, 1, 1};
    static real8 n2[3] = {1, -1, -1};
    static real8 n3[3] = {-1, 1, -1};
    static real8 n4[3] = {-1, -1, 1};

	if (ParallelTest(n1,n)) {
        b1[0]=0;		b1[1]=1;		b1[2]=-1;
        dir1[0]=8.164966e-01;	dir1[1]=-4.082483e-01;	dir1[2]=-4.082483e-01;
        b2[0]=1;		b2[1]=-1;		b2[2]=0;
        dir2[0]=-4.082483e-01;	dir2[1]=-4.082483e-01;	dir2[2]=8.164966e-01;
        b3[0]=-1;		b3[1]=0;		b3[2]=1;
        dir3[0]=-4.082483e-01;	dir3[1]=8.164966e-01;	dir3[2]=-4.082483e-01;
	} else if (ParallelTest(n2,n)) {
        b1[0]=0;		b1[1]=1;		b1[2]=-1;
        dir1[0]=8.164966e-01;	dir1[1]=4.082483e-01;	dir1[2]=4.082483e-01;
        b2[0]=1;		b2[1]=0;		b2[2]=1;
        dir2[0]=-4.082483e-01;	dir2[1]=-8.164966e-01;	dir2[2]=4.082483e-01;
        b3[0]=-1;		b3[1]=-1;		b3[2]=0;
        dir3[0]=-4.082483e-01;	dir3[1]=4.082483e-01;	dir3[2]=-8.164966e-01;
	} else if (ParallelTest(n3,n)) {
        b1[0]=1;		b1[1]=0;		b1[2]=-1;
        dir1[0]=4.082483e-01;	dir1[1]=8.164966e-01;	dir1[2]=4.082483e-01;
        b2[0]=0;		b2[1]=1;		b2[2]=1;
        dir2[0]=-8.164966e-01;	dir2[1]=-4.082483e-01;	dir2[2]=4.082483e-01;
        b3[0]=-1;		b3[1]=-1;		b3[2]=0;
        dir3[0]=4.082483e-01;	dir3[1]=-4.082483e-01;	dir3[2]=-8.164966e-01;
	} else { //must be n4
        b1[0]=-1;		b1[1]=1;		b1[2]=0;
        dir1[0]=4.082483e-01;	dir1[1]=4.082483e-01;	dir1[2]=8.164966e-01;
        b2[0]=1;		b2[1]=0;		b2[2]=1;
        dir2[0]=4.082483e-01;	dir2[1]=-8.164966e-01;	dir2[2]=-4.082483e-01;
        b3[0]=0;		b3[1]=-1;		b3[2]=-1;
        dir3[0]=-8.164966e-01;	dir3[1]=4.082483e-01;	dir3[2]=-4.082483e-01;
	}	

	if (ParallelTest(b1,b)) {
		dir[0]=dir1[0];		dir[1]=dir1[1];	dir[2]=dir1[2];
	} else if (ParallelTest(b2,b)) {
		dir[0]=dir2[0];		dir[1]=dir2[1];	dir[2]=dir2[2];
	} else { //must be b3
		dir[0]=dir3[0];		dir[1]=dir3[1];	dir[2]=dir3[2];
	}

	return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       CrossSlipStress
 *
 *      Description:    Calculate the stress components necessary to calculate
 *						the cross-slip energy barrier.
 *
 *-------------------------------------------------------------------------*/
void CrossSlipStress(Home_t *home, Node_t *node, real8 stress[3][3], 
						real8 linedir[3], real8 burg[3], real8 newplane[3],
						real8 *Ssg, real8 *Seg, real8 *Sscs, real8 *Secs) {

    real8   linesign;
	real8	glideNormal[3], csNormal[3];
    real8   glideTraction[3], csTraction[3];
	real8	glidePosEscaigDir[3], csPosEscaigDir[3];

	Param_t *param;
	param = home->param; 

/*
 *  Vector tangent to dislocation line - the Schmid direction
 *  Use the Burgers vector since the local line direction may not be exactly
 *  screw (due to tolerance on screw detection).
 */
    linesign = (real8) Sign(DotProduct(linedir,burg));
    linedir[0] = linesign*burg[X]; 
    linedir[1] = linesign*burg[Y]; 
    linedir[2] = linesign*burg[Z]; 
	NormalizeVec(linedir);

/*
 *  Grab the glide plane normals that points inside the Thompson tetrahedron. We
 *  need to use this normal because it decides the positive Escaig direction.
 */
    glideNormal[X] = node->nx[0];
    glideNormal[Y] = node->ny[0];
    glideNormal[Z] = node->nz[0];
    ThompsonInwardPlane(glideNormal);
    NormalizeVec(glideNormal);
    csNormal[X] = newplane[X];
	csNormal[Y] = newplane[Y];
	csNormal[Z] = newplane[Z];
    ThompsonInwardPlane(csNormal);
    NormalizeVec(csNormal);

/*
 *  Calculate the traction vectors on the two planes.
 */
    Matrix33Vector3Multiply(stress, glideNormal, glideTraction);
    Matrix33Vector3Multiply(stress, csNormal, csTraction);

/*
 *  Determine the positive Escaig directions.
 */
	FindEscaigDir(glideNormal, burg, glidePosEscaigDir);
	FindEscaigDir(csNormal, burg, csPosEscaigDir);

/*
 *  Determine the relevant traction (stress) components.
 */
    *Seg = DotProduct(glideTraction,glidePosEscaigDir);
    *Ssg = DotProduct(glideTraction,linedir);
 
    *Secs = DotProduct(csTraction,csPosEscaigDir);
    *Sscs = DotProduct(csTraction,linedir);

	return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       EnergyBarrierCrossSlipFCC
 *
 *      Description:    Calculate the energy barrier for cross-slip
 *						in an FCC metal given the Schmid and Escaig
 *						stresses on the glide and cross-slip planes.
 *
 *		This energy barrier function is based on atomistic simulations of Ni. 
 *		See Kuykendall and Cai (2015) for a full explanation.
 *
 *-------------------------------------------------------------------------*/

void EnergyBarrierCrossSlipFCC(Home_t *home, Node_t *node, real8 burg[3], 
								real8 linedir[3], real8 newplane[3], real8 *Eb){

	real8	stress[3][3];
	real8	Sgs, Sge, Scs, Sce;	
	real8	A, T0, p, q, Cge, Ccs, Cce, Tstar, arg; 
	real8	EbFE, EbFl;

	Param_t *param;
	param = home->param;

/*
 *	Calculate stress state at the node and extract key stress components. 
 */
	GetFieldPointStress(home, node->x, node->y, node->z, stress);

	stress[0][0] += param->appliedStress[0];
	stress[1][1] += param->appliedStress[1];
	stress[2][2] += param->appliedStress[2];  
	stress[0][1] += param->appliedStress[5];  
	stress[0][2] += param->appliedStress[4];  
	stress[1][2] += param->appliedStress[3];
	stress[1][0] = stress[0][1];   
	stress[2][0] = stress[0][2];    
	stress[2][1] = stress[1][2];

	CrossSlipStress(home, node, stress, linedir, burg, newplane, &Sgs, &Sge, 
						&Scs, &Sce);

/*
 *	Calculate energy barriers.
 */
	//Freidel-Escaig mechanism
	A = 2.244;			//eV
   	T0 = 5.575e9;		//Pa
   	p = 0.7856;
   	q = 1.937;
   	Cge = -1.679;
   	Ccs = -0.4247;
   	Cce = 0.8823;
	Tstar = Cge*Sge+fabs(Ccs*Scs)+Cce*Sce;
	if (Tstar<0) {
		EbFE = A;
	} else {
		arg = 1-pow(Tstar/T0,p);
		if (arg<0.0) {
			EbFE = 0.0;
		} else {
			EbFE = A*pow(arg,q);
		}	
	}

	//Fleischer mechanism
   	A = 2.428;			//eV
   	T0 = 3.092e9;		//Pa
   	p = 1.270;
   	q = 1.767;
   	Cge = -0.9566;
   	Ccs = -0.7555;
   	Cce = 0.6915;
	Tstar = Cge*Sge+fabs(Ccs*Scs)+Cce*Sce;
	if (Tstar<0) {
		EbFl = A;
	} else {
		arg = 1-pow(Tstar/T0,p);
		if (arg<0.0) {
			EbFl = 0.0;
		} else {
			EbFl = A*pow(arg,q);
		}		
	}

/*
 *	Take the minimum of the energy barrier between the two mechanisms.
 */
	*Eb = MIN(EbFE,EbFl);

	return;
}

