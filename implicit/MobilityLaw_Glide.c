/**************************************************************************
 *
 *  Function    : Mobility_Glide
 *  Author      : Ryan Sills 
 *  Description : Generic Glide Mobility Law 
				  Based on Mobility_FCC_0, allows any plane to be a glide plane
 *                Each line has a glide plane from input
 *                and it never changes
 *                If node flag == 7, node velocity is zero
 *
 *  Returns:  0 on success
 *            1 if velcoity could not be determined
 *
 ***************************************************************************/
#ifdef _IMPLICIT

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif


int Mobility_Glide(Home_t *home, Node_t *node)
{
    int numNonZeroLenSegs = 0;
    Param_t *param;
    real8 VelxNode, VelyNode, VelzNode, Veldotlcr;
    int i, j, nc, nconstraint, nlc;
    real8 normX[100], normY[100], normZ[100], normx[100], normy[100], normz[100];
    real8 a, b;
    real8 dx, dy, dz, lx, ly, lz, lr, LtimesB;
    real8 lcx, lcy, lcz, normdotlc;
    Node_t *nbr;
    real8 MobScrew, MobEdge, Mob;
    real8 bx, by, bz, br, dangle;
    real8 nForce[3];

    param = home->param;

    MobScrew = param->MobScrew;
    MobEdge  = param->MobEdge;
    
    nc = node->numNbrs;
    
/*
 *  If node is 'pinned' in place by constraints, or the node has any arms 
 *  with a burgers vector that has explicitly been set to be sessile (via
 *  control file inputs), the node may not be moved so just zero the velocity
 *  and return
 */
    if ((node->constraint == PINNED_NODE) ||
        NodeHasSessileBurg(home, node))
    {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }

/*
 *  It's possible this function was called for a node which had only zero-
 *  length segments (during SplitSurfaceNodes() for example).  If that is
 *  the case, just set the velocity to zero and return.
 */
	int	zeroLenSegs[100]; //Added by R. B. Sills, 1/9/14, used below
    for (i = 0; i < nc; i++) {
        if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) continue;
        dx = node->x - nbr->x;
        dy = node->y - nbr->y;
        dz = node->z - nbr->z;
		zeroLenSegs[i] = 1;
        if ((dx*dx + dy*dy + dz*dz) > 1.0e-12) {
            numNonZeroLenSegs++;
			zeroLenSegs[i] = 0;
        }
    }

    if (numNonZeroLenSegs == 0) {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }


    /* copy glide plane constraints, but only for segments that have a 
		nonzero length (zero length segments must be ignored) */
    for(i=0;i<nc;i++)
    {
		//if (zeroLenSegs[i]==1) {
			//normX[i] = normx[i] = 0;
		    //normY[i] = normy[i] = 0;
		    //normZ[i] = normz[i] = 0;
		//} else {
		    normX[i] = normx[i] = node->nx[i];
		    normY[i] = normy[i] = node->ny[i];
		    normZ[i] = normz[i] = node->nz[i];
		//}
	
 

/*
 *      If needed, rotate the glide plane normals from the
 *      laboratory frame to the crystal frame.
 */
        if (param->useLabFrame) {
            real8 normTmp[3] = {normX[i], normY[i], normZ[i]};
            real8 normRot[3];

            Matrix33Vector3Multiply(home->rotMatrixInverse, normTmp, normRot);

            normX[i] = normRot[0]; normY[i] = normRot[1]; normZ[i] = normRot[2];
            normx[i] = normRot[0]; normy[i] = normRot[1]; normz[i] = normRot[2];
        }

    }
    
    /* normalize glide plane normal vectors */
    for(i=0;i<nc;i++)
    {
        a=sqrt(normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i]);

        if(a>0)
        {
            normX[i]/=a;
            normY[i]/=a;
            normZ[i]/=a;

            normx[i]/=a;
            normy[i]/=a;
            normz[i]/=a;
        } 
    }


    /* Find independent glide constraints */ 
    nconstraint = nc;
    for(i=0;i<nc;i++)
    {
        for(j=0;j<i;j++)
        {
            Orthogonalize(normX+i,normY+i,normZ+i,normX[j],normY[j],normZ[j]);
        }
        //if((normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i])<FFACTOR_ORTH)
		if((normX[i]*normX[i]+normY[i]*normY[i]+normZ[i]*normZ[i])<1e-4)
        {
            normX[i] = normY[i] = normZ[i] = 0;
            nconstraint--;
        }
    }

    /* find total dislocation length times drag coefficent (LtimesB)*/
    LtimesB=0;
    for(j=0;j<nc;j++)
    {
        if ((nbr=GetNeighborNode(home,node,j)) == (Node_t *)NULL) continue;
        dx=nbr->x - node->x;
        dy=nbr->y - node->y;
        dz=nbr->z - node->z;
        ZImage (param, &dx, &dy, &dz) ;

/*
 *      If needed, rotate the line sense from the laboratory frame to
 *      the crystal frame.
 */
        if (param->useLabFrame) {
            real8 dTmp[3] = {dx, dy, dz};
            real8 dRot[3];

            Matrix33Vector3Multiply(home->rotMatrixInverse, dTmp, dRot);

            dx = dRot[0]; dy = dRot[1]; dz = dRot[2];
        }

        lr=sqrt(dx*dx+dy*dy+dz*dz);
        
        if (lr==0)
        { /* zero arm segment can happen after node split 
           * it is OK to have plane normal vector == 0
           * Skip (do nothing)
           */
        }
        else 
        {
           if((node->nx[j]==0)&&(node->ny[j]==0)&&(node->nz[j]==0))
           {
              printf("Mobility_Glide: (%d,%d) glide plane norm = 0\n"
                     "for segment with nonzero length lr = %e!\n",
                     node->myTag.domainID, node->myTag.index, lr);
           }

           lx=dx/lr; ly=dy/lr; lz=dz/lr;

           bx = node->burgX[j];
           by = node->burgY[j];
           bz = node->burgZ[j];
/*
 *         If needed, rotate the burgers vector from the laboratory frame to
 *         the crystal frame.
 */
           if (param->useLabFrame) {
               real8 bTmp[3] = {bx, by, bz};
               real8 bRot[3];

               Matrix33Vector3Multiply(home->rotMatrixInverse, bTmp, bRot);

               bx = bRot[0]; by = bRot[1]; bz = bRot[2];
           }

           br = sqrt(bx*bx+by*by+bz*bz);
           bx/=br; by/=br; bz/=br; /* unit vector along Burgers vector */

           dangle = fabs(bx*lx+by*ly+bz*lz);
           Mob=MobEdge+(MobScrew-MobEdge)*dangle;

           LtimesB+=(lr/Mob);
	}
    }
    LtimesB/=2;

    nForce[0] = node->fX;
    nForce[1] = node->fY;
    nForce[2] = node->fZ;

/*
 *  If needed, rotate the force vector from the laboratory frame to the
 *  crystal frame
 */
    if (param->useLabFrame) {
        real8 rotForce[3];
        Matrix33Vector3Multiply(home->rotMatrixInverse, nForce, rotForce);
        VECTOR_COPY(nForce, rotForce);
    }

    /* Velocity is simply proportional to total force per unit length */
    VelxNode = nForce[0]/LtimesB;
    VelyNode = nForce[1]/LtimesB;
    VelzNode = nForce[2]/LtimesB;


    /* Orthogonalize with glide plane constraints */
    for(i=0;i<nc;i++)
    {
        if((normX[i]!=0)||(normY[i]!=0)||(normZ[i]!=0))
        {
	     Orthogonalize(&VelxNode,&VelyNode,&VelzNode,
                         normX[i],normY[i],normZ[i]);
        }
    }


/*
 *  If needed, rotate the velocity vector back to the laboratory frame
 *  from the crystal frame
 */
    if (param->useLabFrame) {
        real8 vTmp[3] = {VelxNode, VelyNode, VelzNode};
        real8 vRot[3];

        Matrix33Vector3Multiply(home->rotMatrix, vTmp, vRot);

        VelxNode = vRot[0];
        VelyNode = vRot[1];
        VelzNode = vRot[2];
    }

    node->vX = VelxNode;
    node->vY = VelyNode;
    node->vZ = VelzNode;

    return(0);
}
#endif
