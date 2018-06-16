#include <stdio.h>
#include <math.h>

#ifndef _CYGWIN
#include <complex.h>
#endif
#include "Home.h"

#ifdef _SHDIS
#include "Comm.h"
#include "SH.h"


/*-------------------------------------------------------------------------
 *
 *      Function:     InitSegSigbRem
 *      Description:  Zero out the sigbRem value for all segments.  This
 *                    should only be needed when computing the FEM image
 *                    stress with the FastMultipole(FMM) code enabled (with
 *                    FMM disabled, ComputeSegSigbRem() gets called and
 *                    will handle the needed initialization).
 *
 *-----------------------------------------------------------------------*/
void InitSegSigbRem(Home_t *home, int reqType)
{
        int  i, j;
        Node_t *node;

/*
 *      Initialize sigbRem for all segments attached to native
 *      nodes.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;

/*
 *              If we're only doing a partial force calc, don't update
 *              forces if it is not attached to a node
 *              marked for update.
 */
	    for (j = 0; j < node->numNbrs; j++) {
	      node->sigbRem[3*j  ] = 0.0;
	      node->sigbRem[3*j+1] = 0.0;
	      node->sigbRem[3*j+2] = 0.0;
	    }
	}

/*
 *      Initialize sigbRem for all segments attached to ghost
 *      nodes.
 */
        node = home->ghostNodeQ;

        while (node != (Node_t *)NULL) {

	  for (j = 0; j < node->numNbrs; j++) {
	    node->sigbRem[3*j  ] = 0.0;
	    node->sigbRem[3*j+1] = 0.0;
	    node->sigbRem[3*j+2] = 0.0;
	  }
	  
	  node = node->next;
        }
	
        return;
}


void ComputeSHSegSigbRem(Home_t *home, int reqType)
{
        int     inode, armID1, armID2;
        Node_t  *node, *nbr;
	int     ti, nodeIsOwner;

/*
 *      Loop though all native nodes
 */
        for (inode = 0; inode < home->newNodeKeyPtr; inode++) 
	  {
            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
	    }
	    
/*
 *          loop thru the native node's arms.  If the segment is owned
 *          by the neighboring node, don't do the calculation in this
 *          iteration of the loop.
 */
            for (ti = 0; ti < node->numNbrs; ti++) {
		
                nbr = GetNeighborNode(home, node, ti);

                if ((nodeIsOwner = NodeOwnsSeg(home, node, nbr) == 0)) {
                    continue;
                }

/*
 *              If we're only calculating the sigb for a portion of the
 *              nodes, skip this segment if neither node is flagged for
 *              a force update.
 */
                if ((reqType == PARTIAL) &&
                    (((node->flags & NODE_RESET_FORCES) == 0) &&
                     ((nbr->flags  & NODE_RESET_FORCES) == 0))) {
                    continue;
                }

		armID1 = GetArmID(home, node, nbr);
		armID2 = GetArmID(home, nbr, node);
		
		ComputeSH1SegSigbRem(home, node, 
				     nbr, armID1, armID2);
	      }
	    }
	
        return;

}

int GetNumPhysicNodes(Home_t *home)
{
        int    inode;
        Node_t  *node;
/*
 *      Loop though all native nodes
 */
        int nnodes = 0;
        for (inode = 0; inode < home->newNodeKeyPtr; inode++) 
	      {
            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
	           }
            nnodes++;
        }
        return nnodes;
}

#if 0
// coord = 'x1, x2, x3... \\ y1, y2, y3, ...\\z1, z2, z3, ...'
void GetNodeList(Home_t *home, double* XXX, double* YYY, double* ZZZ)
{
        assert(XXX != NULL); //let python handle memory

        int    i,j, inode;
        Node_t  *node;
/*
 *      Loop though all native nodes
 */
        int nnodes = 0;
        for (inode = 0; inode < home->newNodeKeyPtr; inode++) 
	      {
            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
	           }
            XXX[nnode] = node->x;
            YYY[nnode] = node->y;
            ZZZ[nnode] = node->z;
            nndoes++;
	      }
        return;
}
#else
// coord = 'x1, y1, z1, x2, y2, z2...'
void GetNodeList(Home_t *home, double* XYZ)
{
        assert(XYZ != NULL); //let python handle memory

        int    i,j, inode;
        Node_t  *node;
/*
 *      Loop though all native nodes
 */
        int nnodes = 0;
        for (inode = 0; inode < home->newNodeKeyPtr; inode++) 
	      {
            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
	           }
            j = nnodes * 3;
            XYZ[j]   = node->x;
            XYZ[j+1] = node->y;
            XYZ[j+2] = node->z;
            nnodes++;
	      }
        return;
}
#endif

void SH_calc_stress(Home_t *home, double* stress)
{
  assert(stress != NULL);
        int    i,j, inode;
        Node_t  *node;
#if 0        
        for (inode = 0; inode < home->newNodeKeyPtr; inode++) 
	      {
            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
	           }
            for (i = 0;i<3;i++) for (j = 0;j<3;j++)
                printf("shstress[%d][%d] = %g\n",i,j,stress[inode*9+ i *3 + j ] );
        }
#endif
/*
 *      Loop though all native nodes
 */
        int nnodes = 0;
        for (inode = 0; inode < home->newNodeKeyPtr; inode++) 
	      {
            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
	           }
            for (i = 0;i<3;i++) for (j = 0;j<3;j++)
            {
                //assert(node->SHstress != NULL);
                node->SHstress[i][j] = stress[nnodes *3*3 + i*3 + j];
            }
            nnodes ++;
	      }
        return;
}


void ComputeSH1SegSigbRem(Home_t *home,
			  Node_t *node, Node_t *nbr,
			  int armID1, int armID2)
{

  int     i, j;
  real8   bx, by, bz;
  real8   dx, dy, dz;
  real8   x1, y1, z1;
  real8   sigb1, sigb2, sigb3;
  real8   r[3], totRemSig[3][3];
  
  Param_t *param;
  param = home->param;

  x1 = node->x; 
  y1 = node->y; 
  z1 = node->z;
  
  /*
   *      Get the midpoint of the segment
   */
  dx = nbr->x - x1; 
  dy = nbr->y - y1; 
  dz = nbr->z - z1;
  
  ZImage(param, &dx, &dy, &dz);

  r[0] = x1 + dx*0.5;
  r[1] = y1 + dy*0.5;
  r[2] = z1 + dz*0.5;
  
  /*
   *      Add SH image stress
   */
  
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
    {
//      printf("sh1stress[%d][%d] = %g\n",i,j,node->SHstress[i][j] );
//      printf("sh2stress[%d][%d] = %g\n",i,j,nbr->SHstress[i][j] );
      totRemSig[i][j] = 0.5*(node->SHstress[i][j] + nbr->SHstress[i][j]);
    }
      
  
  /*
   *      Calculate the segment's sig dot b and store it in both nodes
   */
  bx = node->burgX[armID1];
  by = node->burgY[armID1];
  bz = node->burgZ[armID1];
  
  sigb1 = totRemSig[0][0]*bx +
          totRemSig[0][1]*by +
          totRemSig[0][2]*bz;
  sigb2 = totRemSig[1][0]*bx +
          totRemSig[1][1]*by +
          totRemSig[1][2]*bz;
  sigb3 = totRemSig[2][0]*bx +
          totRemSig[2][1]*by +
          totRemSig[2][2]*bz;

  node->sigbRem[3*armID1  ] += sigb1;
  node->sigbRem[3*armID1+1] += sigb2;
  node->sigbRem[3*armID1+2] += sigb3;
  
  nbr->sigbRem[3*armID2  ] += sigb1;
  nbr->sigbRem[3*armID2+1] += sigb2;
  nbr->sigbRem[3*armID2+2] += sigb3;

 
  return;
  
}

#endif
