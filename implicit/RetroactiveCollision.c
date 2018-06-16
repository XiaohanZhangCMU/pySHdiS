/*****************************************************************************
 *
 *      Module:         RetroactiveCollision.c
 *      Description:    This module contains various functions used
 *                      for detecting various types of collisions
 *                      (segment/segment, node/node, zipping) and
 *                      dealing with those collisions.  These functions
 *                      are specific to the type 3 collision handling
 *                      which determines if a collision occurred
 *                      retroactively based on the current and old
 *                      nodal positions. It specifically used the
 *                      interval bisection method to make this
 *                      determination. See Sills and Cai (2014) for
 *                      details.
 *
 *                      Each domain handles local collisions and
 *                      distributes the necessary topological changes
 *                      to remote domains via the same mechanism
 *                      employed by remesh.
 *
 *                      NOTE: Certain strict rules govern what topological
 *                      changes are permitted based on noda and segment
 *                      ownership.  See comments at the beginning of
 *                      the module Topology.c for details of the
 *                      rule.  Additional restrictions may be implemented
 *                      for the collision handling; see code below
 *                      for details of any such restrictions.
 *
 *      Included functions:
 *
 *          FindCollisionPoint()
 *          FindCollisionPointAndTime()
 *          PredictiveCollisions()
 *
 *****************************************************************************/
#ifdef _RETROCOLLISIONS

#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "Home.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"

void SortNodes(Home_t *home, real8 *maxsep);


#define	 MAX_COLS 10000  //Maximum number of collisions allowed in one time step
static real8 COL_TOL = 0.01;  //tolerance for detecting collisions in fraction of time step
static real8 NO_COLLISIONS_TIME = 20.0; //Collision time that denotes no collisions occurred


static int dbgDom;

/*---------------------------------------------------------------------------
 *
 *      Function:       FindCollisionPoint
 *      Description:    This function attempts to select a collision
 *                      point on a plane common to the two nodes.
 *
 *      Arguments:
 *          node1, node2   pointers to the two node structures
 *          x, y, z        pointers to locations in which to return
 *                         the coordinates of the point at which the
 *                         two nodes should be collided.
 *
 *-------------------------------------------------------------------------*/
static void FindCollisionPoint(Home_t *home, Node_t *node1, Node_t *node2,
                               real8 *x, real8 *y, real8 *z)
{
        int     i, j, m, n;
        int     conditionsmet, Nsize;
        int     planeDefined;
        real8   L, invL, tmp;
        real8   norm, invnorm;
        real8   n1mag2, n2mag2, eps;
        real8   dx, dy, dz;
        real8   n1x, n1y, n1z;
        real8   n2x, n2y, n2z;
        real8   dirx, diry, dirz;
        real8   p1[3], p2[3];
        real8   plane[3], vector[3];
        real8   newplanecond, npc2, onemnpc4, detN;
        real8   Nmat[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
        real8   Matrix[6][6], invMatrix[6][6];
        real8   V[6], result[6];
        Node_t  *nbrNode;
        Param_t *param;

        param = home->param;

        eps = 1.0e-12;

        p1[0] = node1->x;
        p1[1] = node1->y;
        p1[2] = node1->z;

        p2[0] = node2->x;
        p2[1] = node2->y;
        p2[2] = node2->z;

        PBCPOSITION(param, p1[0], p1[1], p1[2], &p2[0], &p2[1], &p2[2]);

/*
 *      If a node is a 'fixed' node it can't be relocated, so use
 *      that node's coordinates as the collision point
 */
        if (node1->constraint == PINNED_NODE) {
            *x = p1[0];
            *y = p1[1];
            *z = p1[2];
            return;
        } else if (node2->constraint == PINNED_NODE) {
            *x = p2[0];
            *y = p2[1];
            *z = p2[2];
            return;
        }


       newplanecond = 0.875;
       npc2         = newplanecond * newplanecond;

       tmp          = 1.0 - newplanecond;
       onemnpc4     = tmp * tmp * tmp * tmp;

       vector[0] = 0.0;
       vector[1] = 0.0;
       vector[2] = 0.0;

       Nsize = 0;

       for (i = 0; i < node1->numNbrs; i++) {

           if (Nsize < 3) {

               nbrNode = GetNeighborNode(home, node1, i);

               if (nbrNode == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               dx = p1[0] - nbrNode->x;
               dy = p1[1] - nbrNode->y;
               dz = p1[2] - nbrNode->z;

               ZImage(param, &dx, &dy, &dz);

               L = sqrt(dx*dx + dy*dy + dz*dz);
               invL = 1.0 / L;
               dirx = dx * invL;
               diry = dy * invL;
               dirz = dz * invL;

               xvector(dirx, diry, dirz, node1->burgX[i], node1->burgY[i],
                       node1->burgZ[i], &n1x, &n1y, &n1z);

               xvector(dirx, diry, dirz, node1->vX, node1->vY, node1->vZ,
                       &n2x, &n2y, &n2z);

               n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
               n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

               planeDefined = 0;

               if (n2mag2 > eps) {
/*
 *                 Preference for plane defined by l cross v
 */
                   norm = sqrt(n2mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n2x * invnorm;
                   plane[1] = n2y * invnorm;
                   plane[2] = n2z * invnorm;
                   planeDefined = 1;
               } else if (n1mag2 > eps) {
/*
 *                 Preference for plane defined by l cross b
 */
                   norm = sqrt(n1mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n1x * invnorm;
                   plane[1] = n1y * invnorm;
                   plane[2] = n1z * invnorm;
                   planeDefined = 1;
               } else if (param->enforceGlidePlanes) {
/*
 *                 Use segment's defined glide plane if glide planes enforced
 */
                   plane[0] = node1->nx[i];
                   plane[1] = node1->ny[i];
                   plane[2] = node1->nz[i];
                   NormalizeVec(plane);
                   planeDefined = 1;
               }


               if (planeDefined) {

                   switch (Nsize) {
                   case 0:
                       conditionsmet = 1;
                       break;

                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;

                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];

                      detN = Matrix33Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }

                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, p1);
                       Nsize++;
                   }
               }
           }
       }

       for (i = 0; i < node2->numNbrs; i++) {

           if (Nsize < 3) {

               nbrNode = GetNeighborNode(home, node2, i);

               if (nbrNode == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               dx = p2[0] - nbrNode->x;
               dy = p2[1] - nbrNode->y;
               dz = p2[2] - nbrNode->z;

               ZImage(param, &dx, &dy, &dz);

               L = sqrt(dx*dx + dy*dy + dz*dz);
               invL = 1.0 / L;
               dirx = dx * invL;
               diry = dy * invL;
               dirz = dz * invL;

               xvector(dirx, diry, dirz, node2->burgX[i], node2->burgY[i],
                       node2->burgZ[i], &n1x, &n1y, &n1z);

               xvector(dirx, diry, dirz, node2->vX, node2->vY, node2->vZ,
                       &n2x, &n2y, &n2z);

               n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
               n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

               planeDefined = 0;

               if (n2mag2 > eps) {
/*
 *                 Preference for plane defined by l cross v
 */
                   norm = sqrt(n2mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n2x * invnorm;
                   plane[1] = n2y * invnorm;
                   plane[2] = n2z * invnorm;
                   planeDefined = 1;
               } else if (n1mag2 > eps) {
/*
 *                 Preference for plane defined by l cross b
 */
                   norm = sqrt(n1mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n1x * invnorm;
                   plane[1] = n1y * invnorm;
                   plane[2] = n1z * invnorm;
                   planeDefined = 1;
               } else if (param->enforceGlidePlanes) {
/*
 *                 Use segment's defined glide plane if glide planes enforced
 */
                   plane[0] = node2->nx[i];
                   plane[1] = node2->ny[i];
                   plane[2] = node2->nz[i];
                   NormalizeVec(plane);
                   planeDefined = 1;
               }

               if ((n1mag2 > eps) || (n2mag2 > eps)) {
                   switch (Nsize) {
                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;
                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];
                      detN = Matrix33Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }
 
                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, p2);
                       Nsize++;
                   }
               }
           }
       }

/*
 *         Upper left 3X3 of Matrix is identity matrix.
 *         Matrix rows 3 thru 3+(Nsize-1) colums 0 thru 2 are Nmat.
 *         Matrix columns 3 thru 3+(Nsize-1) rows 0 thru 2 are transpose of Nmat.
 *         All remaining elements are zeroed.
 */
       for (i = 0; i < 6; i++) {
           for (j = 0; j < 6; j++) {
               Matrix[i][j] = 0.0;
           }
       }

       Matrix[0][0] = 1.0;
       Matrix[1][1] = 1.0;
       Matrix[2][2] = 1.0;

       for (i = 0; i < Nsize; i++) {
           for (j = 0; j < 3; j++) {
               Matrix[3+i][j] = Nmat[i][j];
               Matrix[j][3+i] = Nmat[i][j];
           }
       }

       V[0] = 0.5 * (p1[0] + p2[0]);
       V[1] = 0.5 * (p1[1] + p2[1]);
       V[2] = 0.5 * (p1[2] + p2[2]);
       V[3] = vector[0];
       V[4] = vector[1];
       V[5] = vector[2];

       Nsize += 3;

       MatrixInvert((real8 *)Matrix, (real8 *)invMatrix, Nsize, 6);
       MatrixMult((real8 *)invMatrix, Nsize, Nsize, 6,
                  (real8 *)V, 1, 1,
                  (real8 *)result, 1);

       *x = result[0];
       *y = result[1];
       *z = result[2];

       return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       AdjustCollisionPoint
 *      Description:    This function attempts to select a collision
 *                      point on a plane common to the two nodes given a
 *						collision point. This is a modified version of the 
 *						old function FindCollisionPoint. It uses Lagrange 
 *						multipliers to find the point that minimizes the
 *						distance to the cPoint while satisfying all glide
 *						constraints.
 *
 *		Modified by R. B. Sills, 1/7/14
 *
 *      Arguments:
 *			cPoint[3]	   vector containing collision point found by the
						   collision detection algorithm
 *          node1, node2   pointers to the two node structures
 *          x, y, z        pointers to locations in which to return
 *                         the coordinates of the point at which the
 *                         two nodes should be collided.
 *
 *-------------------------------------------------------------------------*/
static void AdjustCollisionPoint(Home_t *home, Node_t *node1, Node_t *node2,
                               real8 cPoint[], real8 *x, real8 *y, real8 *z)
{
        int     i, j, m, n;
        int     conditionsmet, Nsize;
        int     planeDefined;
        real8   L, invL, tmp;
        real8   norm, invnorm;
        real8   n1mag2, n2mag2, eps;
        real8   dx, dy, dz;
        real8   n1x, n1y, n1z;
        real8   n2x, n2y, n2z;
        real8   dirx, diry, dirz;
        real8   p1[3], p2[3];
        real8   plane[3], vector[3];
        real8   newplanecond, npc2, onemnpc4, detN;
        real8   Nmat[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
        real8   Matrix[6][6], invMatrix[6][6];
        real8   V[6], result[6];
        Node_t  *nbrNode;
        Param_t *param;

        param = home->param;

        eps = 1.0e-12;

        p1[0] = node1->x;
        p1[1] = node1->y;
        p1[2] = node1->z;

        p2[0] = node2->x;
        p2[1] = node2->y;
        p2[2] = node2->z;

        PBCPOSITION(param, p1[0], p1[1], p1[2], &p2[0], &p2[1], &p2[2]);
		PBCPOSITION(param, p1[0], p1[1], p1[2], &cPoint[0], &cPoint[1], &cPoint[2]);

/*
 *      If a node is a 'fixed' node it can't be relocated, so use
 *      that node's coordinates as the collision point
 */
        if (node1->constraint == PINNED_NODE) {
            *x = p1[0];
            *y = p1[1];
            *z = p1[2];
            return;
        } else if (node2->constraint == PINNED_NODE) {
            *x = p2[0];
            *y = p2[1];
            *z = p2[2];
            return;
        }


       newplanecond = 0.875;
       npc2         = newplanecond * newplanecond;

       tmp          = 1.0 - newplanecond;
       onemnpc4     = tmp * tmp * tmp * tmp;

       vector[0] = 0.0;
       vector[1] = 0.0;
       vector[2] = 0.0;

/*
 *	   Loop over all arms of both nodes and determine the number of unique glide
 * 	   planes. The Lagrange multipler method can handle at most 3. Populate the
 *     Nmat matrix and vector with this information to be used later.
 */
       Nsize = 0;

       for (i = 0; i < node1->numNbrs; i++) {

           if (Nsize < 3) {

               nbrNode = GetNeighborNode(home, node1, i);

               if (nbrNode == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               dx = p1[0] - nbrNode->x;
               dy = p1[1] - nbrNode->y;
               dz = p1[2] - nbrNode->z;

               ZImage(param, &dx, &dy, &dz);

               L = sqrt(dx*dx + dy*dy + dz*dz);
               invL = 1.0 / L;
               dirx = dx * invL;
               diry = dy * invL;
               dirz = dz * invL;

               xvector(dirx, diry, dirz, node1->burgX[i], node1->burgY[i],
                       node1->burgZ[i], &n1x, &n1y, &n1z);

               xvector(dirx, diry, dirz, node1->vX, node1->vY, node1->vZ,
                       &n2x, &n2y, &n2z);

               n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
               n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

               planeDefined = 0;

               if (param->enforceGlidePlanes) {
/*
 *                 Use segment's defined glide plane if glide planes enforced
 */
                   plane[0] = node1->nx[i];
                   plane[1] = node1->ny[i];
                   plane[2] = node1->nz[i];
                   NormalizeVec(plane);
                   planeDefined = 1;
               }else if (n2mag2 > eps) {
/*
 *                 Preference for plane defined by l cross v
 */
                   norm = sqrt(n2mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n2x * invnorm;
                   plane[1] = n2y * invnorm;
                   plane[2] = n2z * invnorm;
                   planeDefined = 1;
               } else if (n1mag2 > eps) {
/*
 *                 Preference for plane defined by l cross b
 */
                   norm = sqrt(n1mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n1x * invnorm;
                   plane[1] = n1y * invnorm;
                   plane[2] = n1z * invnorm;
                   planeDefined = 1;
               } 

               if (planeDefined) {

                   switch (Nsize) {
                   case 0:
                       conditionsmet = 1;
                       break;

                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;

                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];

                      detN = Matrix33Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }

                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, p1);
                       Nsize++;
                   }
               }
           }
       }

       for (i = 0; i < node2->numNbrs; i++) {

           if (Nsize < 3) {

               nbrNode = GetNeighborNode(home, node2, i);

               if (nbrNode == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               dx = p2[0] - nbrNode->x;
               dy = p2[1] - nbrNode->y;
               dz = p2[2] - nbrNode->z;

               ZImage(param, &dx, &dy, &dz);

               L = sqrt(dx*dx + dy*dy + dz*dz);
               invL = 1.0 / L;
               dirx = dx * invL;
               diry = dy * invL;
               dirz = dz * invL;

               xvector(dirx, diry, dirz, node2->burgX[i], node2->burgY[i],
                       node2->burgZ[i], &n1x, &n1y, &n1z);

               xvector(dirx, diry, dirz, node2->vX, node2->vY, node2->vZ,
                       &n2x, &n2y, &n2z);

               n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
               n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

               planeDefined = 0;

               if (param->enforceGlidePlanes) {
/*
 *                 Use segment's defined glide plane if glide planes enforced
 */
                   plane[0] = node2->nx[i];
                   plane[1] = node2->ny[i];
                   plane[2] = node2->nz[i];
                   NormalizeVec(plane);
                   planeDefined = 1;
               }else if (n2mag2 > eps) {
/*
 *                 Preference for plane defined by l cross v
 */
                   norm = sqrt(n2mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n2x * invnorm;
                   plane[1] = n2y * invnorm;
                   plane[2] = n2z * invnorm;
                   planeDefined = 1;
               } else if (n1mag2 > eps) {
/*
 *                 Preference for plane defined by l cross b
 */
                   norm = sqrt(n1mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n1x * invnorm;
                   plane[1] = n1y * invnorm;
                   plane[2] = n1z * invnorm;
                   planeDefined = 1;
               } 

               if (planeDefined) {
                   switch (Nsize) {
                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;
                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];
                      detN = Matrix33Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }
 
                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, p2);
                       Nsize++;
                   }
               }
           }
       }

/*
 *         Upper left 3X3 of Matrix is identity matrix.
 *         Matrix rows 3 thru 3+(Nsize-1) colums 0 thru 2 are Nmat.
 *         Matrix columns 3 thru 3+(Nsize-1) rows 0 thru 2 are transpose of Nmat.
 *         All remaining elements are zeroed.
 */
       for (i = 0; i < 6; i++) {
           for (j = 0; j < 6; j++) {
               Matrix[i][j] = 0.0;
           }
       }

       Matrix[0][0] = 1.0;
       Matrix[1][1] = 1.0;
       Matrix[2][2] = 1.0;

       for (i = 0; i < Nsize; i++) {
           for (j = 0; j < 3; j++) {
               Matrix[3+i][j] = Nmat[i][j];
               Matrix[j][3+i] = Nmat[i][j];
           }
       }

       V[0] = cPoint[0];
       V[1] = cPoint[1];
       V[2] = cPoint[2];
       V[3] = vector[0];
       V[4] = vector[1];
       V[5] = vector[2];

       Nsize += 3;

       MatrixInvert((real8 *)Matrix, (real8 *)invMatrix, Nsize, 6);
       MatrixMult((real8 *)invMatrix, Nsize, Nsize, 6,
                  (real8 *)V, 1, 1,
                  (real8 *)result, 1);

       *x = result[0];
       *y = result[1];
       *z = result[2];

       return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       AdjustMergePoint
 *      Description:    This function attempts to select a collision
 *                      point on a plane common to all arms of the given 
						node. This is a modified version of the 
 *						old function FindCollisionPoint. It uses Lagrange 
 *						multipliers to find the point that minimizes the
 *						distance to the node's position while satisfying all glide
 *						constraints.
 *
 *		Modified by R. B. Sills, 1/16/14
 *
 *      Arguments:
 *          node		   pointers to the node structure
 *          x, y, z        pointers to locations in which to return
 *                         the coordinates of the point to which the
 *                         node should be moved.
 *
 *-------------------------------------------------------------------------*/
static void AdjustMergePoint(Home_t *home, Node_t *node, real8 *x, 
															real8 *y, real8 *z)
{
        int     i, j, m, n;
        int     conditionsmet, Nsize;
        int     planeDefined;
        real8   L, invL, tmp;
        real8   norm, invnorm;
        real8   n1mag2, n2mag2, eps;
        real8   dx, dy, dz;
        real8   n1x, n1y, n1z;
        real8   n2x, n2y, n2z;
        real8   dirx, diry, dirz;
        real8   p[3], parm[3];
        real8   plane[3], vector[3];
        real8   newplanecond, npc2, onemnpc4, detN;
        real8   Nmat[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
        real8   Matrix[6][6], invMatrix[6][6];
        real8   V[6], result[6];
        Node_t  *nbrNode;
        Param_t *param;

        param = home->param;

        eps = 1.0e-12;

        p[0] = node->x;
        p[1] = node->y;
        p[2] = node->z;

/*
 *      If the node is a 'fixed' node it can't be relocated, so use
 *      that node's coordinates as the collision point
 */
        if (node->constraint == PINNED_NODE) {
            *x = p[0];
            *y = p[1];
            *z = p[2];
            return;
        } 


       newplanecond = 0.875;
       npc2         = newplanecond * newplanecond;

       tmp          = 1.0 - newplanecond;
       onemnpc4     = tmp * tmp * tmp * tmp;

       vector[0] = 0.0;
       vector[1] = 0.0;
       vector[2] = 0.0;

/*
 *	   Loop over all arms of the node and determine the number of unique glide
 * 	   planes. The Lagrange multipler method can handle at most 3. Populate the
 *     Nmat matrix and vector with this information to be used later.
 */
       Nsize = 0;

       for (i = 0; i < node->numNbrs; i++) {

           if (Nsize < 3) {

               nbrNode = GetNeighborNode(home, node, i);

               if (nbrNode == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               dx = p[0] - nbrNode->x;
               dy = p[1] - nbrNode->y;
               dz = p[2] - nbrNode->z;

			   parm[0] = nbrNode->x; parm[1] = nbrNode->y; parm[2] = nbrNode->z;

			   PBCPOSITION(param, p[0], p[1], p[2], &parm[0], &parm[1], &parm[2]);

               ZImage(param, &dx, &dy, &dz);

               L = sqrt(dx*dx + dy*dy + dz*dz);
               invL = 1.0 / L;
               dirx = dx * invL;
               diry = dy * invL;
               dirz = dz * invL;

               xvector(dirx, diry, dirz, node->burgX[i], node->burgY[i],
                       node->burgZ[i], &n1x, &n1y, &n1z);

               xvector(dirx, diry, dirz, node->vX, node->vY, node->vZ,
                       &n2x, &n2y, &n2z);

               n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
               n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

               planeDefined = 0;

               if (param->enforceGlidePlanes) {
/*
 *                 Use segment's defined glide plane if glide planes enforced
 */
                   plane[0] = node->nx[i];
                   plane[1] = node->ny[i];
                   plane[2] = node->nz[i];
                   NormalizeVec(plane);
                   planeDefined = 1;
               }else if (n2mag2 > eps) {
/*
 *                 Preference for plane defined by l cross v
 */
                   norm = sqrt(n2mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n2x * invnorm;
                   plane[1] = n2y * invnorm;
                   plane[2] = n2z * invnorm;
                   planeDefined = 1;
               } else if (n1mag2 > eps) {
/*
 *                 Preference for plane defined by l cross b
 */
                   norm = sqrt(n1mag2);
                   invnorm = 1.0 / norm;
                   plane[0] = n1x * invnorm;
                   plane[1] = n1y * invnorm;
                   plane[2] = n1z * invnorm;
                   planeDefined = 1;
               } 

               if (planeDefined) {

                   switch (Nsize) {
                   case 0:
                       conditionsmet = 1;
                       break;

                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;

                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];

                      detN = Matrix33Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }

                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, parm);
                       Nsize++;
                   }
               }
           }
       }

/*
 *         Upper left 3X3 of Matrix is identity matrix.
 *         Matrix rows 3 thru 3+(Nsize-1) colums 0 thru 2 are Nmat.
 *         Matrix columns 3 thru 3+(Nsize-1) rows 0 thru 2 are transpose of Nmat.
 *         All remaining elements are zeroed.
 */
       for (i = 0; i < 6; i++) {
           for (j = 0; j < 6; j++) {
               Matrix[i][j] = 0.0;
           }
       }

       Matrix[0][0] = 1.0;
       Matrix[1][1] = 1.0;
       Matrix[2][2] = 1.0;

       for (i = 0; i < Nsize; i++) {
           for (j = 0; j < 3; j++) {
               Matrix[3+i][j] = Nmat[i][j];
               Matrix[j][3+i] = Nmat[i][j];
           }
       }

       V[0] = p[0];
       V[1] = p[1];
       V[2] = p[2];
       V[3] = vector[0];
       V[4] = vector[1];
       V[5] = vector[2];

       Nsize += 3;

       MatrixInvert((real8 *)Matrix, (real8 *)invMatrix, Nsize, 6);
       MatrixMult((real8 *)invMatrix, Nsize, Nsize, 6,
                  (real8 *)V, 1, 1,
                  (real8 *)result, 1);

       *x = result[0];
       *y = result[1];
       *z = result[2];

       return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       IntervalCollision
 *      Description:    Apply the interval halving method.                      
 *
 *      Arguments:
 *          cPoint      Will be set to the coordinates of the collision
 *                      point if the two segments will collide.
 *          cTime       Set to the time (in units of current timestep
 *                      length (deltaTT)) at which the two segments will
 *                      collide, or negative if not colliding.
 *          minDist2    Min distance between segments, squared.
 *          L1, L2		Normalized distance to minimum distance points
 *			x/y/zclosei	Coordinates of the ith minimum distance point
 *			ddist2dt	Time rate of change of the square of the minimum
 *						separation distance
 *
 *-------------------------------------------------------------------------*/
static void IntervalCollision(Home_t *home, real8 *L1out, real8 *L2out, real8 *hitTime,
									  real8 x1old, real8 y1old, real8 z1old,
                                      real8 x2old, real8 y2old, real8 z2old,
                                      real8 x3old, real8 y3old, real8 z3old,
                                      real8 x4old, real8 y4old, real8 z4old,
                                      real8 x1, real8 y1, real8 z1,
                                      real8 x2, real8 y2, real8 z2,
                                      real8 x3, real8 y3, real8 z3,
                                      real8 x4, real8 y4, real8 z4, 
									  real8 start, real8 finish, real8 drA,
									  real8 drB, real8 *Dout, int loc, 
									  real8 *Dmin, real8 *minTime, 
									  real8 *L1Min, real8 *L2Min, int *n)

{

		real8 	eps = 1.0e-10;
		real8	tol = 1.0e-2;
		real8	x10, y10, z10, x20, y20, z20, x30, y30, z30, x40, y40, z40; 
		real8	x1f, y1f, z1f, x2f, y2f, z2f, x3f, y3f, z3f, x4f, y4f, z4f;
		real8	Df, D0, Df2, D02;
		real8	L10, L20, L1f, L2f;
		real8	ddist2dt;
		real8	hitTimePass, DminPass, minTimePass, L1MinPass, L2MinPass;
		real8	dint, intdrA, intdrB, maxDistance, mid;
		int		nmax = 200;

/*
 *		Increment the iteration counter.
 */
		(*n)++;
		if (*n>nmax) {
			*L1out = 0.0;
			*L2out = 0.0;
			*Dout = 0.0;
			*hitTime = NO_COLLISIONS_TIME;//;
			return;
		}

/*		
 *		Move the nodes to the positions at the beginning and end of the interval.
 */
		x10 = x1old+start*(x1-x1old);
		y10 = y1old+start*(y1-y1old);
		z10 = z1old+start*(z1-z1old);
		x20 = x2old+start*(x2-x2old);
		y20 = y2old+start*(y2-y2old);
		z20 = z2old+start*(z2-z2old);
		x30 = x3old+start*(x3-x3old);
		y30 = y3old+start*(y3-y3old);
		z30 = z3old+start*(z3-z3old);
		x40 = x4old+start*(x4-x4old);
		y40 = y4old+start*(y4-y4old);
		z40 = z4old+start*(z4-z4old);

		x1f = x1old+finish*(x1-x1old);
		y1f = y1old+finish*(y1-y1old);
		z1f = z1old+finish*(z1-z1old);
		x2f = x2old+finish*(x2-x2old);
		y2f = y2old+finish*(y2-y2old);
		z2f = z2old+finish*(z2-z2old);
		x3f = x3old+finish*(x3-x3old);
		y3f = y3old+finish*(y3-y3old);
		z3f = z3old+finish*(z3-z3old);
		x4f = x4old+finish*(x4-x4old);
		y4f = y4old+finish*(y4-y4old);
		z4f = z4old+finish*(z4-z4old);

/*		
 *		Calculate min separation distance at start and finish.  
 *		Depending on which end of the interval was passed through, we
 *		only need to calculate one. 
 */
		if (loc==1) {
			Df = *Dout;
			L1f = *L1out;	L2f = *L2out;
			GetMinDist2(x10, y10, z10, 0.0, 0.0, 0.0,
				        x20, y20, z20, 0.0, 0.0, 0.0,
				        x30, y30, z30, 0.0, 0.0, 0.0,
				        x40, y40, z40, 0.0, 0.0, 0.0,
						&D02, &ddist2dt, &L10, &L20);
			D0 = sqrt(D02);
		} else {
			D0 = *Dout;
			L10 = *L1out;	L20 = *L2out;
			GetMinDist2(x1f, y1f, z1f, 0.0, 0.0, 0.0,
				        x2f, y2f, z2f, 0.0, 0.0, 0.0,
				        x3f, y3f, z3f, 0.0, 0.0, 0.0,
				        x4f, y4f, z4f, 0.0, 0.0, 0.0,
						&Df2, &ddist2dt, &L1f, &L2f);
			Df = sqrt(Df2); 
		}

/*		
 *		Update the minimum distance, locations, and time overall all of time
 */
		*Dmin = MIN(*Dmin, D0);
		*Dmin = MIN(*Dmin, Df);
		if (*Dmin==D0) {
			*minTime = start;
			*L1Min = L10;
			*L2Min = L20;
		} else if (*Dmin==Df) {
			*minTime = finish;
			*L1Min = L1f;
			*L2Min = L2f;
		}			

/*		
 *		Calculate max distance moved by both segments over the interval. 
 */
		dint = finish-start;
		intdrA = dint*drA;
		intdrB = dint*drB;
		maxDistance = intdrA+intdrB;
/*
 *		If the max distance is less than the minimum distance at the beginning
 *		or end of the interval, there will be no collision in this interval.
 */
		if (maxDistance<D0 || maxDistance<Df) {
			*L1out = 0.0;
			*L2out = 0.0;
			*Dout = 0.0;
			*hitTime = NO_COLLISIONS_TIME;
			return;
		}

/*
 *		If the interval is smaller than the tolerance, consider the segments
 *		as collided.  Take the average of the minimun distance positions
 *		and time interval for output.
 */
		if (finish-start < COL_TOL) {
			*L1out = (L10+L1f)/2.0;
			*L2out = (L20+L2f)/2.0;
			*Dout = 0.0;
			*hitTime = (start+finish)/2.0;
			return;
		}

/*		
 *		Bisect the interval and continue the recursion.
 */
		mid = (start + finish)/2.0;

/*		
 *		Check for a collision in the first half of the inteveral first.
 *		If no collision there, then look in the second half.
 */

		IntervalCollision(home, &L10, &L20, hitTime,
							x1old, y1old, z1old, x2old, y2old, z2old,
	                        x3old, y3old, z3old, x4old, y4old, z4old,
	                        x1, y1, z1, x2, y2, z2,
	                        x3, y3, z3, x4, y4, z4, 
							start, mid, drA, drB, &D0, 0,
							Dmin, minTime, L1Min, L2Min, n);

		if (*n>nmax) {
			*L1out = 0.0;
			*L2out = 0.0;
			*Dout = 0.0;
			*hitTime = NO_COLLISIONS_TIME;
			return;
		}

		if (*hitTime <=1) {
			*L1out = L10;
			*L2out = L20;
			*Dout = D0;
			return;
		}

		IntervalCollision(home, &L1f, &L2f, hitTime,
							x1old, y1old, z1old, x2old, y2old, z2old,
	                        x3old, y3old, z3old, x4old, y4old, z4old,
	                        x1, y1, z1, x2, y2, z2,
	                        x3, y3, z3, x4, y4, z4, 
							mid, finish, drA, drB, &Df, 1,
							Dmin, minTime, L1Min, L2Min, n);

		*L1out = L1f;
		*L2out = L2f;
		*Dout = Df;

		return; 
}


/*---------------------------------------------------------------------------
 *
 *      Function:       DetectCollisionBisection
 *      Description:    Given the nodal positions from the last and current
 *						time step, determine if two segments collided by 
 *						comparing the directions of the vectors connecting
 *						their closest points and their minimum separation
 *						distances.                     
 *
 *      Arguments:
 *          cPoint      Will be set to the coordinates of the collision
 *                      point if the two segments will collide.
 *          cTime       Set to the time (in units of current timestep
 *                      length (deltaTT)) at which the two segments will
 *                      collide, or negative if not colliding.
 *          minDist2    Min distance between segments, squared.
 *          L1, L2		Normalized distance to minimum distance points
 *			x/y/zclosei	Coordinates of the ith minimum distance point
 *			ddist2dt	Time rate of change of the square of the minimum
 *						separation distance
 *
 *-------------------------------------------------------------------------*/
static void DetectCollisionBisection(Home_t *home, real8 cPoint[3], int *collided,
									  real8 *L1, real8 *L2, real8 *cTime,
									  real8 x1old, real8 y1old, real8 z1old,
                                      real8 x2old, real8 y2old, real8 z2old,
                                      real8 x3old, real8 y3old, real8 z3old,
                                      real8 x4old, real8 y4old, real8 z4old,
                                      real8 x1, real8 y1, real8 z1,
                                      real8 x2, real8 y2, real8 z2,
                                      real8 x3, real8 y3, real8 z3,
                                      real8 x4, real8 y4, real8 z4, 
									  real8 *dist2)                                  
{
		real8	drx, dry, drz, dr1, dr2, dr3, dr4, drA, drB;
		real8	hitTime, L1c, L2c, dist2old, distold; 
		real8	x1c, y1c, z1c, x2c, y2c, z2c, x3c, y3c, z3c, x4c, y4c, z4c; 
		real8	cPoint1[3], cPoint2[3];
		real8	distMin, minTime, L1Min, L2Min, distinit;
		int		n;
		Param_t	*param;
		
		param = home->param; 

/*
 *		Calculate the maximum amount moved by either segment.
 */
		drx = x1-x1old;	dry = y1-y1old;	drz = z1-z1old;
		dr1 = sqrt(drx*drx+dry*dry+drz*drz);

		drx = x2-x2old;	dry = y2-y2old;	drz = z2-z2old;
		dr2 = sqrt(drx*drx+dry*dry+drz*drz);

		drA = MAX(dr1,dr2);

		drx = x3-x3old;	dry = y3-y3old;	drz = z3-z3old;
		dr3 = sqrt(drx*drx+dry*dry+drz*drz);

		drx = x4-x4old;	dry = y4-y4old;	drz = z4-z4old;
		dr4 = sqrt(drx*drx+dry*dry+drz*drz);

		drB = MAX(dr3,dr4);		

/*		
 *		Package the minimum separation distance and positions in the old 
 *		configuration for passing to the collision function.
 */
		distold = sqrt(*dist2);
        distinit = distold;
		L1c = *L1;
		L2c = *L2;
		distMin = distold;
		minTime = 0.0;
		L1Min = *L1;
		L2Min = *L2;
		n = 0;

/*
 *		Begin the interval halving method...
 */
		IntervalCollision(home, &L1c, &L2c, &hitTime, 
							x1old, y1old, z1old,
		                    x2old, y2old, z2old,
		                    x3old, y3old, z3old,
		                    x4old, y4old, z4old,
		                    x1, y1, z1,
		                    x2, y2, z2,
		                    x3, y3, z3,
		                    x4, y4, z4, 0.0, 1.0, drA, drB, 
							&distold, 0, 
							&distMin, &minTime, &L1Min, &L2Min, &n);

/*
 *      If the collision time is less than the collision tolerance
 *      and the segments started at a separation distance greater than
 *      the annihilation distance, it is likely a collision did not
 *      actually occur. Prevent the collision from being detected,
 *      unless distMin<rann, as tested below.
 */
        if (hitTime<COL_TOL && distinit>param->rann) {
            hitTime = NO_COLLISIONS_TIME;
        }

/*
 *		If the returned hit time (as a fraction of the last time step)
 *		is less than or equal to 1, there was a collision.  Otherwise, if the minimum
 * 		distance found while looking for a collision was less than the annihilation
 *		radius, collide the segments at that point and time.
 */
		if (hitTime<=1) {  
			x1c = x1old+hitTime*(x1-x1old);
			y1c = y1old+hitTime*(y1-y1old);
			z1c = z1old+hitTime*(z1-z1old);
			x2c = x2old+hitTime*(x2-x2old);
			y2c = y2old+hitTime*(y2-y2old);
			z2c = z2old+hitTime*(z2-z2old);
			x3c = x3old+hitTime*(x3-x3old);
			y3c = y3old+hitTime*(y3-y3old);
			z3c = z3old+hitTime*(z3-z3old);
			x4c = x4old+hitTime*(x4-x4old);
			y4c = y4old+hitTime*(y4-y4old);
			z4c = z4old+hitTime*(z4-z4old);

			cPoint1[0] = x1c+L1c*(x2c-x1c);
			cPoint1[1] = y1c+L1c*(y2c-y1c);
			cPoint1[2] = z1c+L1c*(z2c-z1c);
			cPoint2[0] = x3c+L2c*(x4c-x3c);
			cPoint2[1] = y3c+L2c*(y4c-y3c);
			cPoint2[2] = z3c+L2c*(z4c-z3c);
			cPoint[0] = (cPoint1[0]+cPoint2[0])/2.0;
			cPoint[1] = (cPoint1[1]+cPoint2[1])/2.0;
			cPoint[2] = (cPoint1[2]+cPoint2[2])/2.0;
			*collided = 1;
			*cTime = hitTime;
			*L1 = L1c;
			*L2 = L2c;
		} else if (distMin<param->rann) {
			x1c = x1old+minTime*(x1-x1old);
			y1c = y1old+minTime*(y1-y1old);
			z1c = z1old+minTime*(z1-z1old);
			x2c = x2old+minTime*(x2-x2old);
			y2c = y2old+minTime*(y2-y2old);
			z2c = z2old+minTime*(z2-z2old);
			x3c = x3old+minTime*(x3-x3old);
			y3c = y3old+minTime*(y3-y3old);
			z3c = z3old+minTime*(z3-z3old);
			x4c = x4old+minTime*(x4-x4old);
			y4c = y4old+minTime*(y4-y4old);
			z4c = z4old+minTime*(z4-z4old);

			cPoint1[0] = x1c+L1Min*(x2c-x1c);
			cPoint1[1] = y1c+L1Min*(y2c-y1c);
			cPoint1[2] = z1c+L1Min*(z2c-z1c);
			cPoint2[0] = x3c+L2Min*(x4c-x3c);
			cPoint2[1] = y3c+L2Min*(y4c-y3c);
			cPoint2[2] = z3c+L2Min*(z4c-z3c);
			cPoint[0] = (cPoint1[0]+cPoint2[0])/2.0;
			cPoint[1] = (cPoint1[1]+cPoint2[1])/2.0;
			cPoint[2] = (cPoint1[2]+cPoint2[2])/2.0;
			*collided = 1;
			*cTime = minTime;
			*L1 = L1Min;
			*L2 = L2Min;
		} else {
			*collided = 0;
			*cTime = NO_COLLISIONS_TIME;
			*L1 = 0.0;
			*L2 = 0.0;
			cPoint[0] = 0.0;	cPoint[1] = 0.0;	cPoint[2] = 0.0;
		}
		return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       TestGlidePlanes
 *      Description:    Given a node, tests all of its arms for a glide  plane
 *						violation, and returns a 0 or a 1 through the passTest
 *						variable. Glide plane is considered violated if the arm
 *						deviates from the plane by more than the length 
 *						tol_test*rann.
 *
 *                      Since a node may already be violating glide planes, and
 *                      we don't want to bias these nodes against collisions,
 *                      we only test whether the pre-existing violation is made
 *                      worse.                    
 *
 *      Arguments:
 *          *home		home structure pointer
 *			*node		pointer for node to test
 *			*passtest	integer pointer for variable that returns test status
 *
 *-------------------------------------------------------------------------*/	

void TestGlidePlanes(Home_t *home, Node_t *node, real8 p[3], int *passtest)
{
		int		j;
		real8	nx, ny, nz, nbrx, nbry, nbrz, drx, dry, drz, drlen;
		real8	tol_test, tol_len, dottest, violen_old, violen_max;
		Node_t	*nbr;
		Param_t	*param;

		param = home->param;
		tol_test = 1e-2;  //tolerance for plane test
		tol_len = 1e-10;  //tolerance for "zero-length" segment


        violen_max = tol_test*param->rann;
		for (j=0; j < node->numNbrs; j++) {
			nx = node->nx[j];
			ny = node->ny[j];
			nz = node->nz[j];
			Normalize(&nx,&ny,&nz);
			nbr = GetNeighborNode(home, node, j);
			nbrx = nbr->x;
			nbry = nbr->y;
			nbrz = nbr->z;

            //First test the current position.
			PBCPOSITION(param, node->x, node->y, node->z, &nbrx, &nbry, &nbrz);
			drx = node->x-nbrx;
			dry = node->y-nbry;
			drz = node->z-nbrz;

            drlen = sqrt(drx*drx+dry*dry+drz*drz);
			if (drlen<tol_len) continue;

			Normalize(&drx,&dry,&drz);
			dottest = fabs(nx*drx+ny*dry+nz*drz);
            violen_old = dottest*drlen;

            //Then test the new position.
			PBCPOSITION(param, p[0], p[1], p[2], &nbrx, &nbry, &nbrz);
			drx = p[0]-nbrx;
			dry = p[1]-nbry;
			drz = p[2]-nbrz;

			drlen = sqrt(drx*drx+dry*dry+drz*drz);
			if (drlen<tol_len) continue;

			Normalize(&drx,&dry,&drz);
			dottest = fabs(nx*drx+ny*dry+nz*drz);

            //If the violation increases by more than violen_max,
            //it fails the test.
			if ((dottest*drlen - violen_old) > violen_max) {
				*passtest = 0;
				return;					
			}
		}
		*passtest = 1;

		return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    Linked list structures for collisions
 *
 *-----------------------------------------------------------------------*/
typedef struct _collisionitem {
		Node_t  *node1, *node2, *node3, *node4;
		real8   cTime;
		real8   cPoint[3];
		real8   L1old, L2old;
		real8   p1[3], p2[3], p3[3], p4[3];
		int     id, done;
		struct _collisionitem *next;
} CollisionItem;

typedef struct _collisionlist {
		int size;
		CollisionItem *head;
} CollisionList;

typedef struct _pairitem {
		Node_t *nodes[4];
		struct _pairitem *next;
} PairItem;

typedef struct _pairlist {
		int size;
		PairItem *head;
} PairList;

/*------------------------------------------------------------------------
 *
 *      Function:    Sort4
 *      Description: Sort an array of four integers in ascending order
 *
 *-----------------------------------------------------------------------*/
static __inline__ int SortNodes4(Node_t **n) {
		int i, j;
		for (i = 1; i < 4; i++) {
			Node_t *tmp = n[i];
			for (j = i; j >= 1 && OrderNodes(tmp, n[j-1]) < 0; j--)
				n[j] = n[j-1];
			n[j] = tmp;
		}
                return 0;
}

/*------------------------------------------------------------------------
 *
 *      Function:    ComparePairItems
 *      Description: Compare two pair items
 *
 *-----------------------------------------------------------------------*/
int ComparePairItems(const PairItem *item1, const PairItem *item2)
{
		int i;
		
		for (i = 0; i < 4; i++) {
			if (OrderNodes(item1->nodes[i], item2->nodes[i]) != 0) {
				return OrderNodes(item1->nodes[i], item2->nodes[i]);
			}
		}
		
		return 0;
}

/*------------------------------------------------------------------------
 *
 *      Function:    	FindInsertPairList
 *      Description: 	Find if a pair item already exists in the list and
 *                      insert it in ascending order if this is not the case
 *
 *-----------------------------------------------------------------------*/
int FindInsertPairList(PairList *pairlist, Node_t **nodes)
{
		int       order;
		PairItem *item, *cur;
		
		item = (PairItem*)malloc(sizeof(PairItem));
		item->nodes[0] = nodes[0];
		item->nodes[1] = nodes[1];
		item->nodes[2] = nodes[2];
		item->nodes[3] = nodes[3];
		item->next = NULL;
		
		if (pairlist->head == NULL) {
			pairlist->head = item;
			pairlist->size++;
			return 1;
		}
			
		order = ComparePairItems(item, pairlist->head);
		if (order < 0) {
			item->next = pairlist->head;
			pairlist->head = item;
			pairlist->size++;
			return 1;
		} else if (order == 0) {
			free(item);
			return 0;
		}
		
		cur = pairlist->head;
		while (cur != NULL) {
			
			if (cur->next == NULL) {
				cur->next = item;
				pairlist->size++;
				return 1;
			}
			
			order = ComparePairItems(item, cur->next);
			if (order < 0) {
				item->next = cur->next;
				cur->next = item;
				pairlist->size++;
				return 1;
			} else if (order == 0) {
				free(item);
				return 0;
			}
			
			cur = cur->next;
		}
    return 0;
}

/*------------------------------------------------------------------------
 *
 *      Function:    	InsertCollisionList
 *      Description: 	Insert a new collision item in the list after
 * 						checking if the item already existed.
 *
 *-----------------------------------------------------------------------*/
int InsertCollisionList(CollisionList *colList, PairList *pairlist, 
                        Node_t *node1, Node_t *node2, Node_t *node3, Node_t *node4, 
                        real8 cTime, real8 cPoint[3], real8 L1old, real8 L2old, 
                        real8 p1[3], real8 p2[3], real8 p3[3], real8 p4[3])
{
	
/*		
 * 		Check if a similar collision pair already
 * 		exists in the list and insert it if needed
 */
		Node_t *nodes[4];
		nodes[0] = node1;
		nodes[1] = node2;
		nodes[2] = node3;
		nodes[3] = node4;
		SortNodes4(nodes);
		
		if (!FindInsertPairList(pairlist, nodes)) {
			return 0;
		}

/*		
 * 		Insert the collision pair to the list
 */		
		int i;
		CollisionItem *item, *cur;
		
		item = (CollisionItem*)malloc(sizeof(CollisionItem));
		item->id = colList->size;
		item->node1 = node1;
		item->node2 = node2;
		item->node3 = node3;
		item->node4 = node4;
		item->cTime = cTime;
		item->L1old = L1old;
		item->L2old = L2old;
		for (i = 0; i < 3; i++) {
			item->cPoint[i] = cPoint[i];
			item->p1[i] = p1[i];
			item->p2[i] = p2[i];
			item->p3[i] = p3[i];
			item->p4[i] = p4[i];
		}
		item->done = 0;
		item->next = NULL;
		
		if (colList->head == NULL) {
			colList->head = item;
			colList->size++;
			return 1;
		}
		
		cur = colList->head;
		while (cur->next != NULL) {
			cur = cur->next;
		}
		cur->next = item;
		colList->size++;
		
		return 1;
}

/*------------------------------------------------------------------------
 *
 *      Function:    	UpdateCollisionList
 *      Description: 	Update the collision list by removing every collision
 * 						item that is flagged as done
 *
 *-----------------------------------------------------------------------*/
void UpdateCollisionList(CollisionList *colList)
{
		CollisionItem *head, *cur, *next;
		
		cur = colList->head;
		if (cur == NULL) return;
		
		while (cur->next != NULL) {
			if (cur->next->done == 1) {
				next = cur->next;
				cur->next = next->next;
				free(next);
			} else {
				cur = cur->next;
			}
		}
		
		cur = colList->head;
		if (cur->done == 1) {
			colList->head = cur->next;
			free(cur);
		}
		
		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:		EmptyCollisionList
 *      Description:	Empty the collision list
 *
 *-----------------------------------------------------------------------*/
void EmptyCollisionList(CollisionList *colList)
{
		CollisionItem *cur, *next;
		
		cur = colList->head;
		while (cur != NULL) {
			next = cur->next;
			free(cur);
			cur = next;
		}
		
		colList->head = NULL;
		colList->size = 0;
}

/*------------------------------------------------------------------------
 *
 *      Function:		EmptyPairList
 *      Description:	Empty the collision pair list
 *
 *-----------------------------------------------------------------------*/
void EmptyPairList(PairList *pairlist)
{
		PairItem *cur, *next;
		
		cur = pairlist->head;
		while (cur != NULL) {
			next = cur->next;
			free(cur);
			cur = next;
		}
		
		pairlist->head = NULL;
		pairlist->size = 0;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       RetroactiveCollisions
 *      Description:    Loop though all native nodes, identify segments
 *                      or nodes that should be collided by
 *                      the current domain and handle all such collisions.
 *
 *                      We detect collisions using the nodal positions from the
 *                      previous (store in node->olderxyz) and current time 
 *                      steps using the interval halving method See Sills et al.
 *                      (2015) in MSMSE for details.
 *
 *                      Note the following restrictions on collisions:
 *
 *                      - A 'pinned' node be deleted during a collision
 *
 *      All nodes with arms owned by another domain have previously
 *      been marked as non-deletable nodes for this cycle.
 * 
 *      A linked list is used to efficienctly store, retreive and remove
 *      the collision items as the collision algorithm proceeds.
 *
 *-------------------------------------------------------------------------*/
void RetroactiveCollisions(Home_t *home)
{
        int     i, j, k, q, arm12, arm21, arm34, arm43;
        int     thisDomain, splitStatus, mergeStatus;
        int     armAB, armBA;
        int     globalOp = 1;
        int     close2node1, close2node2, close2node3, close2node4;
        int     splitSeg1, splitSeg2;
        int     cell2Index, nbrCell2Index, nextIndex;
        int     cell2X, cell2Y, cell2Z, cx, cy, cz;
        int     localCollisionCnt, globalCollisionCnt;
        int     collisionConditionIsMet, adjustCollisionPoint;
        real8   dist2, ddist2dt, L1, L2, eps, tol, half, mindist;
		real8	dist2old, ddist2dtold, L1old, L2old;
        real8   cTime, cPoint[3];
        real8   x1, y1, z1, vx1, vy1, vz1;
        real8   x2, y2, z2, vx2, vy2, vz2;
        real8   x3, y3, z3, vx3, vy3, vz3;
        real8   x4, y4, z4, vx4, vy4, vz4;
        real8   x1old, y1old, z1old, vx1old, vy1old, vz1old;
        real8   x2old, y2old, z2old, vx2old, vy2old, vz2old;
        real8   x3old, y3old, z3old, vx3old, vy3old, vz3old;
        real8   x4old, y4old, z4old, vx4old, vy4old, vz4old;
        real8   newx, newy, newz;
        real8   newvx, newvy, newvz;
        real8   pnew[3], burg1[3], burg2[3];
        real8   oldfp0s1[3], oldfp1s1[3], oldfp0s2[3], oldfp1s2[3];
        real8   f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];
        real8   nodeVel[3], newNodeVel[3];
        real8   newPos[3], newVel[3];
        real8   vec1[3], vec2[3];
        real8   vec1old[3], vec2old[3];
        real8   p1[3], p2[3], p3[3], p4[3];
        real8   v1[3], v2[3], v3[3], v4[3];
		real8 	xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3, xc4, yc4, zc4;
		real8 	xclose1, yclose1, zclose1, xclose2, yclose2, zclose2;
		real8 	xc1old, yc1old, zc1old, xc2old, yc2old, zc2old;
		real8 	xc1ish, yc1ish, zc1ish, xc2ish, yc2ish, zc2ish;
		real8	earliestcTime, tmp;
		real8	merge1pos[3], merge1vel[3], merge2pos[3], merge2vel[3];
        real8   mindist2;
		int		morecollisions, colnum, col, alreadydone;
		int		narms1, narms2, passtest, arm;
        Tag_t   oldTag1, oldTag2;
        Node_t  *node1, *node2, *node3, *node4, *tmpNbr;
        Node_t  *mergenode1, *mergenode2, *targetNode;
        Node_t  *splitNode1, *splitNode2;
		Node_t	*tmpNode;
		int 	donesize, ndone;
        Param_t *param;
#ifdef _FEM
        int     resetSurfaceProperties, femSurface[2];
        real8   femSurfaceNorm[3];
#endif

        TimerStart(home, COLLISION_HANDLING);

/*
 *      Allocate memory for collision arrays.
 */
        real8	**done_points = malloc(MAX_COLS*sizeof(real8 *));
		
		for(i=0;i<MAX_COLS;i++){
			done_points[i] = malloc(3 * sizeof(real8));
        }
		
		CollisionList *colList;
		colList = (CollisionList*)malloc(sizeof(CollisionList));
		colList->size = 0;
		colList->head = NULL;
		
		CollisionItem *cur;
		
		PairList *pairlist;
		pairlist = (PairList*)malloc(sizeof(PairList));
        

        thisDomain = home->myDomain;
        param      = home->param;

        eps      = 1.0e-12;
		tol		 = 1e-3;
        half     = 0.5;
        mindist2 = param->rann * param->rann;

        localCollisionCnt = 0;
        globalCollisionCnt = 0;

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom = -1;
#endif

//		Initialize collision variables
		earliestcTime = NO_COLLISIONS_TIME;
		cTime = 10.0;
		donesize = 0;
		ndone = 0;
		colnum = 0;
		morecollisions = 1;
		int passcollisions = -1;

/*
 *		Loop until all collisions are handled.
 */
		while (morecollisions) {
			passcollisions++;
			pairlist->size = 0;
			pairlist->head = NULL;

/*
 *      Start looping through native nodes looking for segments to collide...
 *      We only consider segment pairs that don't share a common node (all
 *      collisions can be handled in this way). 
 */
	    for (i = 0; i < home->newNodeKeyPtr; i++) {

	        if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
	        if (node1->flags & NO_COLLISIONS) continue;
	        
/*
 *          Loop through all cell2s neighboring the node.  Only
 *          nodes in these neighboring cell2s are candidates for
 *          collisions.
 *
 *          NOTE: All nodes are assigned cell2 membership prior to
 *          entering collision handling.  However, nodes added during
 *          node separation and collision handling are not given
 *          cell2 membership for this cycle.  In most cases this
 *          is not a problem since those nodes have usually been
 *          exempted from subsequent collisions this cycle.  There
 *          are certain circumstances, though, in which the exemptions
 *          have not been set.  If we encounter one such node that
 *          has no cell2 membership, skip it for this cycle, or bad
 *          things happen.
 */
	        cell2Index = node1->cell2Idx;
	        if (cell2Index < 0) {
	            continue;
	        }
	        
	        x1 = node1->x;  x1old = node1->olderx;
	        y1 = node1->y;  y1old = node1->oldery;
	        z1 = node1->z;  z1old = node1->olderz;
	        
	        vx1 = node1->vX;  vx1old = node1->oldervX;
	        vy1 = node1->vY;  vy1old = node1->oldervY;
	        vz1 = node1->vZ;  vz1old = node1->oldervZ;
	        
	        PBCPOSITION(param, x1, y1, z1, &x1old, &y1old, &z1old);
	        
	        p1[X] = x1;  v1[X] = vx1;
	        p1[Y] = y1;  v1[Y] = vy1;
	        p1[Z] = z1;  v1[Z] = vz1;

	        
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
	                if (node3->flags & NO_COLLISIONS) continue;

	                if (CollisionNodeOrder(home, &node1->myTag,
	                                       &node3->myTag) >= 0) {
	                    continue;
	                }
					
					if ((node1->myTag.domainID==node3->myTag.domainID)&&
						(node1->myTag.index   ==node3->myTag.index)) {
						continue;
					}
					
					x3 = node3->x;  x3old = node3->olderx;
					y3 = node3->y;  y3old = node3->oldery;
					z3 = node3->z;  z3old = node3->olderz;
					
					vx3 = node3->vX;  vx3old = node3->oldervX;
					vy3 = node3->vY;  vy3old = node3->oldervY;
					vz3 = node3->vZ;  vz3old = node3->oldervZ;
					
					PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
					PBCPOSITION(param, x3, y3, z3, &x3old, &y3old, &z3old);
					
					p3[X] = x3;  v3[X] = vx3;
					p3[Y] = y3;  v3[Y] = vy3;
					p3[Z] = z3;  v3[Z] = vz3;


/*
 *                  Loop over all arms of node1.  Skip any arms that
 *                  terminate at node3 (those hinge arms will be dealt
 *                  with later)
 */
	                for (arm12 = 0; arm12 < node1->numNbrs; arm12++) {

	                    node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);

	                    if (node2 == (Node_t *)NULL) continue;
	                    if (node2->flags & NO_COLLISIONS) continue;


	                    if ((node2->myTag.domainID == node3->myTag.domainID) &&
	                        (node2->myTag.index    == node3->myTag.index   )) {
							continue;
	                    }

#ifdef PARALLEL
						if (CollisionNodeOrder(home, &node1->myTag,
                                               &node2->myTag) > 0) {
                            continue;
                        }

/*
 *                      Segment node1/node2 may only be used in a collision
 *                      if the segment is owned by the current domain.
 */
	                    if (!DomainOwnsSeg(home, OPCLASS_COLLISION,
	                                       thisDomain, &node2->myTag)) {
	                        continue;
	                    }
#endif

						x2 = node2->x;  x2old = node2->olderx;
						y2 = node2->y;  y2old = node2->oldery;
						z2 = node2->z;  z2old = node2->olderz;
						
						vx2 = node2->vX;  vx2old = node2->oldervX;
						vy2 = node2->vY;  vy2old = node2->oldervY;
						vz2 = node2->vZ;  vz2old = node2->oldervZ;
						
						PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
						PBCPOSITION(param, x1old, y1old, z1old, &x2old, &y2old, &z2old);
						
						p2[X] = x2;  v2[X] = vx2;
						p2[Y] = y2;  v2[Y] = vy2;
						p2[Z] = z2;  v2[Z] = vz2;
						
						
/*
 * 						It is possible to have a zero-length segment
 * 						(created by a previous collision).  If we find
 * 						such a segment, do not try to use it in any
 * 						subsequent collisions.
 */
						vec1[X] = x2 - x1;
						vec1[Y] = y2 - y1;
						vec1[Z] = z2 - z1;

						if (DotProduct(vec1, vec1) < 1.0e-20) {
							continue;
						}

						vec1old[X] = x2old - x1old;
						vec1old[Y] = y2old - y1old;
						vec1old[Z] = z2old - z1old;

						if (DotProduct(vec1old, vec1old) < 1.0e-20) {
							continue;
						}

/*
 *                      Loop over all arms of node3.  Skip any segments that
 *                      terminate at node2 (those hinge arms will be dealt
 *                      with later), and skip any segments not owned by
 *                      node3 -- we'll deal with those segments when we
 *                      hit the node owning the segment
 *                      
 */
	                    for (arm34 = 0; arm34 < node3->numNbrs; arm34++) {

/*
 *							Protect against overstepping the size of our arrays.
 */
							if (donesize>=MAX_COLS || colnum>=MAX_COLS) {
								break;
							}

	                        node4 = GetNodeFromTag(home, node3->nbrTag[arm34]);

	                        if (node4 == (Node_t *)NULL) continue;
	                        if (node4->flags & NO_COLLISIONS) continue;

	                        if ((node4->myTag.domainID==node2->myTag.domainID)&&
	                            (node4->myTag.index   ==node2->myTag.index)) {
								continue;
	                        }

	                        if ((node1->myTag.domainID==node4->myTag.domainID)&&
	                            (node1->myTag.index   ==node4->myTag.index)) {
								continue;
	                        }

#ifdef PARALLEL
							if (CollisionNodeOrder(home, &node3->myTag,
                                                   &node4->myTag) > 0) {
                                continue;
                            }

/*
 *                          At this point, segment node3/node4 is owned by
 *                          node3.  If node3 is not native to this domain,
 *                          the segment may not be used in a collision since
 *                          the domain doing to collision must own both
 *                          segments.
 */
	                        if (node3->myTag.domainID != thisDomain) {
	                            continue;
	                        }
							if (!DomainOwnsSeg(home, OPCLASS_COLLISION,
	                                       thisDomain, &node4->myTag)) {
								continue;
							}
							
#endif
							
							collisionConditionIsMet = 0;
							
		                    x4 = node4->x;  x4old = node4->olderx;
		                    y4 = node4->y;  y4old = node4->oldery;
		                    z4 = node4->z;  z4old = node4->olderz;
		                    
		                    vx4 = node4->vX;  vx4old = node4->oldervX;
		                    vy4 = node4->vY;  vy4old = node4->oldervY;
		                    vz4 = node4->vZ;  vz4old = node4->oldervZ;
		                    
		                    PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);
		                    PBCPOSITION(param, x3old, y3old, z3old, &x4old, &y4old, &z4old);
		                    
		                    p4[X] = x4;  v4[X] = vx4;
		                    p4[Y] = y4;  v4[Y] = vy4;
		                    p4[Z] = z4;  v4[Z] = vz4;

/*
 *                          It is possible to have a zero-length segment
 *                          (created by a previous collision).  If we find
 *                          such a segment, do not try to use it in any
 *                          subsequent collisions.
 */
		                    vec2[X] = x4 - x3;
		                    vec2[Y] = y4 - y3;
		                    vec2[Z] = z4 - z3;

		                    if (DotProduct(vec2, vec2) < 1.0e-20) {
		                        continue;
		                    }

		                    vec2old[X] = x4old - x3old;
		                    vec2old[Y] = y4old - y3old;
		                    vec2old[Z] = z4old - z3old;

		                    if (DotProduct(vec2old, vec2old) < 1.0e-20) {
		                        continue;
		                    }

		                    GetMinDist2(x1old, y1old, z1old, vx1old, vy1old, vz1old,
		                               x2old, y2old, z2old, vx2old, vy2old, vz2old,
		                               x3old, y3old, z3old, vx3old, vy3old, vz3old,
		                               x4old, y4old, z4old, vx4old, vy4old, vz4old,
		                               &dist2old, &ddist2dtold, &L1old, &L2old);

							DetectCollisionBisection(home, cPoint, &collisionConditionIsMet,
							                         &L1old, &L2old, &cTime,
							                         x1old, y1old, z1old,
							                         x2old, y2old, z2old,
							                         x3old, y3old, z3old,
							                         x4old, y4old, z4old,
							                         x1, y1, z1, x2, y2, z2,
							                         x3, y3, z3, x4, y4, z4, &dist2old);

/*
 *							If we have already handled a collision with this
 *							collision point, don't add this to the list.
 */
							real8 dx, dy ,dz;
							alreadydone = 0;
							for (j=0; j<ndone; j++) {
								dx = cPoint[X]-done_points[j][X];
								dy = cPoint[Y]-done_points[j][Y];
								dz = cPoint[Z]-done_points[j][Z];
								ZImage(param,&dx,&dy,&dz);
								if ((dx*dx+dy*dy+dz*dz)<0.01*mindist2) {
									alreadydone = 1;
									break;
								} 
							}
							if (alreadydone) continue;
					

/*
 *							If we have a collision, add it to the collision data
 *							arrays.
 */
							if (colnum>=MAX_COLS) break;
							if (collisionConditionIsMet) {
								int insert = InsertCollisionList(colList, pairlist, node1, node2, node3, node4, cTime, 
								                                 cPoint, L1old, L2old, p1, p2, p3, p4);
								if (insert) colnum++;
							}
							
	                    }  /* Loop over node3 arms */
	                }  /* Loop over node1 arms */
	            }  /* while(nextIndex >= 0) */
	        }  /* Loop over neighboring cell2s */
	       }  /* Loop over neighboring cell2s */
	      }  /* Loop over neighboring cell2s */
		} /* Loop over each native node */

/*
 *		Find the earliest collision and handle it. This also tests whether or 
 *		not we have more collisions to handle (potentially).
 */
		earliestcTime = NO_COLLISIONS_TIME;
		morecollisions = 0;
		
		CollisionItem *colitem = NULL;
		cur = colList->head;
		while (cur != NULL) {
			if (cur->cTime < earliestcTime && cur->done == 0) {
				earliestcTime = cur->cTime;
				colitem = cur;
				morecollisions = 1;
			}
			cur = cur->next;
		}

/*
 *		Protect against overstepping the size of our arrays.
 */
		if (donesize>=MAX_COLS || colnum>=MAX_COLS) {
			earliestcTime = NO_COLLISIONS_TIME;
			morecollisions = 0;
			printf("Warning: Collision detection was aborted due to max number"
					"of collisions (MAX_COLS) being reached in RetroactiveCollisions.\n");
		}

/*
 *		Reset the donesize to zero, since we will be doing a new search in the
 *		next iteration.
 */
		donesize = 0;
		EmptyPairList(pairlist);

/*		Unpack the earliest collision data if there was a collision */
		if (earliestcTime<10.0) {

/*			
 *			Unpack earliest collision info.
 */

			node1 = colitem->node1;
			node2 = colitem->node2;
			node3 = colitem->node3;
			node4 = colitem->node4;
			cPoint[X] = colitem->cPoint[X];
			cPoint[Y] = colitem->cPoint[Y];
			cPoint[Z] = colitem->cPoint[Z];
			L1old = colitem->L1old;
			L2old = colitem->L2old;
			x1 = colitem->p1[X];
			y1 = colitem->p1[Y];
			z1 = colitem->p1[Z];	
			x2 = colitem->p2[X];
			y2 = colitem->p2[Y];
			z2 = colitem->p2[Z];
			x3 = colitem->p3[X];
			y3 = colitem->p3[Y];
			z3 = colitem->p3[Z];
			x4 = colitem->p4[X];
			y4 = colitem->p4[Y];
			z4 = colitem->p4[Z];
			
			colitem->done = 1;

            vx1 = node1->vX; vy1 = node1->vY; vz1 = node1->vZ;
            vx2 = node2->vX; vy2 = node2->vY; vz2 = node2->vZ;
            vx3 = node3->vX; vy3 = node3->vY; vz3 = node3->vZ;
            vx4 = node4->vX; vy4 = node4->vY; vz4 = node4->vZ;
		
			arm12 = GetArmID(home,node1,node2);
			arm34 = GetArmID(home,node3,node4);

			done_points[ndone][X] = cPoint[X];
			done_points[ndone][Y] = cPoint[Y];
			done_points[ndone][Z] = cPoint[Z];
			ndone++;
			if (ndone>=MAX_COLS) break;

#ifdef DEBUG_LOG_COLLISIONS
			printf("earlestcTime=%e\n",earliestcTime);				
			printf("Collided at x=%e y=%e z=%e\n",cPoint[X],cPoint[Y],cPoint[Z]);
			printf("Node 1 is %i at x=%e y=%e z=%e\n",node1->myTag.index,x1,y1,z1);
			printf("Node 2 is %i at x=%e y=%e z=%e\n",node2->myTag.index,x2,y2,z2);
			printf("Node 3 is %i at x=%e y=%e z=%e\n",node3->myTag.index,x3,y3,z3);
			printf("Node 4 is %i at x=%e y=%e z=%e\n",node4->myTag.index,x4,y4,z4);
			printf("L1old=%e L2old=%e\n",L1old,L2old); 
#endif

/*
 *			Now before we make any topological changes, we need to do some 
 *			bookkeeping. First, we flag all nodes for NO_COLLISIONS so we don't
 *			redo the whole search on the next iteration. We also need to go 
 *			through the collision data arrays and "delete" any collisions that 
 *			involve any of the nodes in these cell2s, as those collisions may
 *			not be correct (once we handle this collision).
 *
 *			Note that the nodes involved in this particular collision will 
 *			also be unflagged, so that the collision can proceed.		
 */	
	    	for (i = 0; i < home->newNodeKeyPtr; i++) {
	        	if ((tmpNode = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				tmpNode->flags |= NO_COLLISIONS;
			}

			node1->flags &= ~NO_COLLISIONS;
			for (i=0; i<node1->numNbrs; i++) {
				tmpNode = GetNeighborNode(home,node1,i);
				tmpNode->flags &= ~NO_COLLISIONS;
			}
			node2->flags &= ~NO_COLLISIONS;
			for (i=0; i<node2->numNbrs; i++) {
				tmpNode = GetNeighborNode(home,node2,i);
				tmpNode->flags &= ~NO_COLLISIONS;
			}
			node3->flags &= ~NO_COLLISIONS;
			for (i=0; i<node3->numNbrs; i++) {
				tmpNode = GetNeighborNode(home,node3,i);
				tmpNode->flags &= ~NO_COLLISIONS;
			}
			node4->flags &= ~NO_COLLISIONS;
			for (i=0; i<node4->numNbrs; i++) {
				tmpNode = GetNeighborNode(home,node4,i);
				tmpNode->flags &= ~NO_COLLISIONS;
			}

			cur = colList->head;
			while (cur != NULL) {
				if (cur->done == 0) {
					if ((cur->node1->flags & NO_COLLISIONS) == 0 ||
					    (cur->node2->flags & NO_COLLISIONS) == 0 ||
					    (cur->node3->flags & NO_COLLISIONS) == 0 ||
					    (cur->node4->flags & NO_COLLISIONS) == 0) {
						cur->node1->flags &= ~NO_COLLISIONS;
						cur->node2->flags &= ~NO_COLLISIONS;
						cur->node3->flags &= ~NO_COLLISIONS;
						cur->node4->flags &= ~NO_COLLISIONS;
						cur->done = 1;
					}
				}
				cur = cur->next;
			}
/*
 * 			Remove all collisions flagged as done from the list
 */
			UpdateCollisionList(colList);
			
			
/*			Take care of the earliest collision left */

/*
*           Segments are unconnected and colliding.
*           Identify the first node to be merged.  If the
*           collision point is close to one of the nodal
*           endpoints, use that node, otherwise insert a
*           new node in the segment.
*
*           NOTE: The current domain owns node1 but may
*           not own node2.  If it does not own node2, we 
*           cannot allow the collision to use node2
*           even if the collision point is close to 
*           that node.
*/

/*
 *			Determine if either node on segment 1 was close to the collision
 *			point at the time of the collision. If it was within the minimum
 *			segment distance, we use it for the collision. Otherwise, we add
 *			a new node.
 *			
 *			Note that in the past rann was used (instead of minSeg), but this
 *			led to very rightly angled nodes when rann was very small, causing
 *			bad things to happen...
 */
			x1old = node1->olderx; y1old = node1->oldery; z1old = node1->olderz;
			PBCPOSITION(param, x1, y1, z1, &x1old, &y1old, &z1old);

			xc1 = x1old + earliestcTime*(x1-x1old);
			yc1 = y1old + earliestcTime*(y1-y1old);
			zc1 = z1old + earliestcTime*(z1-z1old);

			x2old = node2->olderx; y2old = node2->oldery; z2old = node2->olderz;
			PBCPOSITION(param, x2, y2, z2, &x2old, &y2old, &z2old);

			xc2 = x2old + earliestcTime*(x2-x2old);
			yc2 = y2old + earliestcTime*(y2-y2old);
			zc2 = z2old + earliestcTime*(z2-z2old);

            vec1[X] = cPoint[X] - xc1;
            vec1[Y] = cPoint[Y] - yc1;
            vec1[Z] = cPoint[Z] - zc1;

            vec2[X] = cPoint[X] - xc2;
            vec2[Y] = cPoint[Y] - yc2;
            vec2[Z] = cPoint[Z] - zc2;

            //close2node1 = (DotProduct(vec1,vec1)<mindist2);
            //close2node2 = (DotProduct(vec2,vec2)<mindist2);

            close2node1 = (DotProduct(vec1,vec1)<param->minSeg*param->minSeg);
            close2node2 = (DotProduct(vec2,vec2)<param->minSeg*param->minSeg);

            if ((node2->myTag.domainID != thisDomain) &&
                close2node2) {
				node2->flags |= NO_COLLISIONS;
                continue;
            }

            if (close2node1) {
                mergenode1 = node1;
                splitSeg1 = 0;
            } else if (close2node2) {
                mergenode1 = node2;
                splitSeg1 = 0;
            } else {
                splitSeg1 = 1;
            }

/*
 *          If we need to add a new node to the first segment, do it now.
 */
            if (splitSeg1) {
                 real8 pos0[3], pos1[3];

                 newx = x1 * (1.0-L1old) + x2*L1old;
                 newy = y1 * (1.0-L1old) + y2*L1old;
                 newz = z1 * (1.0-L1old) + z2*L1old;

                 newvx = vx1 * (1.0-L1old) + vx2*L1old;
                 newvy = vy1 * (1.0-L1old) + vy2*L1old;
                 newvz = vz1 * (1.0-L1old) + vz2*L1old;

/*
 *               Estimate resulting forces on all segments involved in the split.
 */
                 arm21 = GetArmID(home, node2, node1);

                 oldfp0s1[X] = node1->armfx[arm12];
                 oldfp0s1[Y] = node1->armfy[arm12];
                 oldfp0s1[Z] = node1->armfz[arm12];

                 oldfp1s1[X] = node2->armfx[arm21];
                 oldfp1s1[Y] = node2->armfy[arm21];
                 oldfp1s1[Z] = node2->armfz[arm21];

                 pos0[X] = x1;   pos1[X] = x2;
                 pos0[Y] = y1;   pos1[Y] = y2;
                 pos0[Z] = z1;   pos1[Z] = z2;

                 pnew[X] = newx;
                 pnew[Y] = newy;
                 pnew[Z] = newz;

                 newNodeVel[X] = newvx;
                 newNodeVel[Y] = newvy;
                 newNodeVel[Z] = newvz;

                 nodeVel[X] = vx1;
                 nodeVel[Y] = vy1;
                 nodeVel[Z] = vz1;

                 burg1[X] = node1->burgX[arm12];
                 burg1[Y] = node1->burgY[arm12];
                 burg1[Z] = node1->burgZ[arm12];

                 FindSubFSeg(home, pos0, pos1, burg1, oldfp0s1,
                             oldfp1s1, pnew, f0seg1, f1seg1,
                             f0seg2, f1seg2);

                 oldTag1 = node1->myTag;
                 oldTag2 = node2->myTag;

                 FoldBox(param, &pnew[X], &pnew[Y], &pnew[Z]);

                 splitStatus = SplitNode(home,
                                         OPCLASS_COLLISION,
                                         node1, pos0, pnew,
                                         nodeVel,
                                         newNodeVel, 1,
                                         &arm12, globalOp,
                                         &splitNode1,
                                         &splitNode2, 0);


/*
 *               If we were unable to split the node
 *               go back to looking for more collision
 *               candidates.
 */
                 if (splitStatus == SPLIT_FAILED) {
                     node1->flags |= NO_COLLISIONS;
                     continue;
                 }

/*
 *				 It is possible that the older positions got messed up
 *               during the split. Also, the new node needs to have its older
 *               position set. Set the older positions here.
 */
				 splitNode1->olderx = x1old;
				 splitNode1->oldery = y1old;
				 splitNode1->olderz = z1old;
				 FoldBox(param,&splitNode1->olderx,&splitNode1->oldery,&splitNode1->olderz);

				 splitNode2->olderx = x1old * (1.0-L1old) + x2old*L1old;
				 splitNode2->oldery = y1old * (1.0-L1old) + y2old*L1old;
				 splitNode2->olderz = z1old * (1.0-L1old) + z2old*L1old;
				 FoldBox(param,&splitNode2->olderx,&splitNode2->oldery,&splitNode2->olderz);

/*
 *               The force estimates above are good enough
 *               for the remainder of this timestep, but
 *               mark the force and velocity data for
 *               some nodes as obsolete so that more
 *               accurate forces will be recalculated
 *               either at the end of this timestep, or
 *               the beginning of the next.
 */
                 mergenode1 = splitNode2;

                 MarkNodeForceObsolete(home, splitNode2);

                 for (q = 0; q < splitNode2->numNbrs; q++) {
                     tmpNbr = GetNodeFromTag(home, splitNode2->nbrTag[q]);
                     if (tmpNbr != (Node_t *)NULL) {
                         tmpNbr->flags |= NODE_RESET_FORCES;
                     }
                 }

/*
 *               Reset nodal forces on nodes involved in the
 *               split.
 */
                 ResetSegForces(home, splitNode1,
                                &splitNode2->myTag,
                                f0seg1[X], f0seg1[Y],
                                f0seg1[Z], 1);

                 ResetSegForces(home, splitNode2,
                                &splitNode1->myTag,
                                f1seg1[X], f1seg1[Y],
                                f1seg1[Z], 1);

                 ResetSegForces(home, splitNode2,
                                &node2->myTag,
                                f0seg2[X], f0seg2[Y],
                                f0seg2[Z], 1);

                 ResetSegForces(home, node2,
                                &splitNode2->myTag,
                                f1seg2[X], f1seg2[Y],
                                f1seg2[Z], 1);

                 (void)EvaluateMobility(home, splitNode1);
                 (void)EvaluateMobility(home, splitNode2);
                 (void)EvaluateMobility(home, node2);

/*
 *               When debugging, dump some info on
 *               topological changes taking place and
 *               the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                if ((dbgDom < 0)||(dbgDom == home->myDomain)) {
                    printf("  Split-1(SegCollision): "
                           "(%d,%d)--(%d,%d) ==> "
                           "(%d,%d)--(%d,%d)--(%d,%d)\n",
                           oldTag1.domainID, oldTag1.index,
                           oldTag2.domainID, oldTag2.index,
                           splitNode1->myTag.domainID,
                           splitNode1->myTag.index,
                           splitNode2->myTag.domainID,
                           splitNode2->myTag.index,
                           node2->myTag.domainID,
                           node2->myTag.index);
                    PrintNode(splitNode1);
                    PrintNode(splitNode2);
                    PrintNode(node2);
                 }
#endif
            }  /* if (splitSeg1) */

/*
 *          Identify the second node to be merged
 *
 *          Note: The current domain owns node3 but may not
 *          own node4.  If it does not own node4, we 
 *          cannot allow the collision to use node4
 *          even if the collision point is close to 
 *          that node.
 */
			x3old = node3->olderx; y3old = node3->oldery; z3old = node3->olderz;
			PBCPOSITION(param, x3, y3, z3, &x3old, &y3old, &z3old);

			xc3 = x3old + earliestcTime*(x3-x3old);
			yc3 = y3old + earliestcTime*(y3-y3old);
			zc3 = z3old + earliestcTime*(z3-z3old);

			x4old = node4->olderx; y4old = node4->oldery; z4old = node4->olderz;
			PBCPOSITION(param, x4, y4, z4, &x4old, &y4old, &z4old);

			xc4 = x4old + earliestcTime*(x4-x4old);
			yc4 = y4old + earliestcTime*(y4-y4old);
			zc4 = z4old + earliestcTime*(z4-z4old);

            vec1[X] = cPoint[X] - xc3;
            vec1[Y] = cPoint[Y] - yc3;
            vec1[Z] = cPoint[Z] - zc3;

            vec2[X] = cPoint[X] - xc4;
            vec2[Y] = cPoint[Y] - yc4;
            vec2[Z] = cPoint[Z] - zc4;

            //close2node3 = (DotProduct(vec1,vec1)<mindist2);
            //close2node4 = (DotProduct(vec2,vec2)<mindist2);

            close2node3 = (DotProduct(vec1,vec1)<param->minSeg*param->minSeg);
            close2node4 = (DotProduct(vec2,vec2)<param->minSeg*param->minSeg);

            if ((node4->myTag.domainID != thisDomain) &&
                close2node4) {
				node4->flags |= NO_COLLISIONS;
                continue;
            }

            if (close2node3) {
                mergenode2 = node3;
                splitSeg2 = 0;
            } else if (close2node4) {
                mergenode2 = node4;
                splitSeg2 = 0;
            } else {
                splitSeg2 = 1;
            }

/*
 *          If we need to add a new node to the second
 *          segment, do it now.
 */
            if (splitSeg2) {
                 real8 pos0[3], pos1[3];

                 newx = x3 * (1.0-L2old) + x4*L2old;
                 newy = y3 * (1.0-L2old) + y4*L2old;
                 newz = z3 * (1.0-L2old) + z4*L2old;

                 newvx = vx3 * (1.0-L2old) + vx4*L2old;
                 newvy = vy3 * (1.0-L2old) + vy4*L2old;
                 newvz = vz3 * (1.0-L2old) + vz4*L2old;

/*
 *               Estimate resulting forces on all segments
 *               involved in the split.
 */
                 arm43 = GetArmID(home, node4, node3);

                 burg1[X] = node3->burgX[arm34];
                 burg1[Y] = node3->burgY[arm34];
                 burg1[Z] = node3->burgZ[arm34];

                 oldfp0s2[X] = node3->armfx[arm34];
                 oldfp0s2[Y] = node3->armfy[arm34];
                 oldfp0s2[Z] = node3->armfz[arm34];

                 oldfp1s2[X] = node4->armfx[arm43];
                 oldfp1s2[Y] = node4->armfy[arm43];
                 oldfp1s2[Z] = node4->armfz[arm43];

                 pos0[X] = x3;   pos1[X] = x4;
                 pos0[Y] = y3;   pos1[Y] = y4;
                 pos0[Z] = z3;   pos1[Z] = z4;

                 pnew[X] = newx;
                 pnew[Y] = newy;
                 pnew[Z] = newz;

                 newNodeVel[X] = newvx;
                 newNodeVel[Y] = newvy;
                 newNodeVel[Z] = newvz;

                 nodeVel[X] = vx3;
                 nodeVel[Y] = vy3;
                 nodeVel[Z] = vz3;

                 FindSubFSeg(home, pos0, pos1, burg1, oldfp0s2,
                             oldfp1s2, pnew, f0seg1, f1seg1,
                             f0seg2, f1seg2);

                 oldTag1 = node3->myTag;
                 oldTag2 = node4->myTag;

                 FoldBox(param, &pnew[X], &pnew[Y], &pnew[Z]);

                 splitStatus = SplitNode(home,
                                        OPCLASS_COLLISION,
                                        node3, pos0,
                                        pnew, nodeVel,
                                        newNodeVel, 1,
                                        &arm34, globalOp,
                                        &splitNode1,
                                        &splitNode2, 0);

/*
 *               If we were unable to split the node
 *               go back to looking for more collision
 *               candidates.
 */
                 if (splitStatus == SPLIT_FAILED) {
                     node3->flags |= NO_COLLISIONS;
                     continue;
                 }

/*
 *				 It is possible that the older positions got messed up
 *               during the split. Also, the new node needs to have its older
 *               position set. Set the older positions here.
 */
				 splitNode1->olderx = x3old;
				 splitNode1->oldery = y3old;
				 splitNode1->olderz = z3old;
				 FoldBox(param,&splitNode1->olderx,&splitNode1->oldery,&splitNode1->olderz);

				 splitNode2->olderx = x3old * (1.0-L2old) + x4old*L2old;
				 splitNode2->oldery = y3old * (1.0-L2old) + y4old*L2old;
				 splitNode2->olderz = z3old * (1.0-L2old) + z4old*L2old;
				 FoldBox(param,&splitNode2->olderx,&splitNode2->oldery,&splitNode2->olderz);

/*
 *               The force estimates above are good enough
 *               for the remainder of this timestep, but
 *               mark the force and velocity data for some
 *               nodes as obsolete so that more accurate
 *               forces will be recalculated either at the
 *               end of this timestep, or the beginning of
 *               the next.
 */
                 mergenode2 = splitNode2;

                 MarkNodeForceObsolete(home, splitNode2);

                 for (q = 0; q < splitNode2->numNbrs; q++) {
                     tmpNbr = GetNodeFromTag(home, splitNode2->nbrTag[q]);
                     if (tmpNbr != (Node_t *)NULL) {
                         tmpNbr->flags |= NODE_RESET_FORCES;
                     }
                 }

/*
 *               Reset nodal forces on nodes involved in the split.
 */
                 ResetSegForces(home, splitNode1,
                                &splitNode2->myTag,
                                f0seg1[X], f0seg1[Y],
                                f0seg1[Z], 1);

                 ResetSegForces(home, splitNode2,
                                &splitNode1->myTag,
                                f1seg1[X], f1seg1[Y],
                                f1seg1[Z], 1);

                 ResetSegForces(home, splitNode2,
                                &node4->myTag,
                                f0seg2[X], f0seg2[Y],
                                f0seg2[Z], 1);

                 ResetSegForces(home, node4,
                                &splitNode2->myTag,
                                f1seg2[X], f1seg2[Y],
                                f1seg2[Z], 1);

                 (void)EvaluateMobility(home, splitNode1);
                 (void)EvaluateMobility(home, splitNode2);
                 (void)EvaluateMobility(home, node4);

/*
 *               When debugging, dump some info on
 *               topological changes taking place and
 *               the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                 if ((dbgDom < 0) ||
                     (dbgDom == home->myDomain)) {
                     printf("  Split-2(SegCollision): "
                           "(%d,%d)--(%d,%d) ==> "
                           "(%d,%d)--(%d,%d)--(%d,%d)\n",
                           oldTag1.domainID,
                           oldTag1.index,
                           oldTag2.domainID,
                           oldTag2.index,
                           splitNode1->myTag.domainID,
                           splitNode1->myTag.index,
                           splitNode2->myTag.domainID,
                           splitNode2->myTag.index,
                           node4->myTag.domainID,
                           node4->myTag.index);
                     PrintNode(splitNode1);
                     PrintNode(splitNode2);
                     PrintNode(node4);
                }
#endif
            }  /* if (splitSeg2) */
			
/*
 *			If both nodes are pinned, we need to unpin one to handle
 *			the collisions - this is needed when precipitates are being
 *			simulated.
 */			
			if (mergenode1->constraint == PINNED_NODE &&
				mergenode2->constraint == PINNED_NODE) {
					mergenode1->constraint &= ~PINNED_NODE;
			}

/*
 *          Adjust the collision point so that it satisfies all glide 
 *			constraints of the involved nodes.			
 */
			AdjustCollisionPoint(home, mergenode1, mergenode2,
                               cPoint, &newPos[X], &newPos[Y], &newPos[Z]);

/*
 *          However, if it looks like the node will be thrown a significant 
 *          distance during the collision (due to glide constraints), 
 *          don't do the collision unless one of the nodes is pinned (in
 *			which case the point was moved to that node).
 */

            vec1[X] = newPos[X] - cPoint[X];
            vec1[Y] = newPos[Y] - cPoint[Y];
            vec1[Z] = newPos[Z] - cPoint[Z];

			if (DotProduct(vec1, vec1) > 16.0 * mindist2 &&
				(mergenode1->constraint != PINNED_NODE &&
				 mergenode2->constraint != PINNED_NODE)) {
                mergenode1->flags |= NO_COLLISIONS;
                mergenode2->flags |= NO_COLLISIONS;
                continue;
            }

/*
 *          Test for glide plane violations with the new collision
 *          point. If we have any, don't handle the collision. 
 *
 *          BAD: This is slightly inconsistent with the idea of 
 *          an annihilation radius, but there's no way around it
 *          it if we are strictly enforcing glide planes...
 */
			int passtest1, passtest2;
			TestGlidePlanes(home, mergenode1, newPos, &passtest1);
			TestGlidePlanes(home, mergenode2, newPos, &passtest2);
			if (passtest1==0 || passtest2==0) {
//printf("Collision aborted due to glide violation!\n");
                mergenode1->flags |= NO_COLLISIONS;
                mergenode2->flags |= NO_COLLISIONS;
                continue;
           }


			FoldBox(param, &newPos[X],&newPos[Y],&newPos[Z]);


#ifdef _FEM
/*
 *         If colliding 2 surface nodes, we may have to
 *         adjust the collision point so it too is on the
 *         surface.
 */
           resetSurfaceProperties = 0;

           if ((mergenode1->constraint == SURFACE_NODE) &&
               (mergenode2->constraint == SURFACE_NODE)) {
               Node_t *seg1Node2, *seg2Node2;

               seg1Node2 = (mergenode1 == node1) ?
                           node2 : node1;
               seg2Node2 = (mergenode2 == node3) ?
                           node4 : node3;

               FEM_AdjustCollisionPoint(mergenode1,
                                        seg1Node2,
                                        mergenode2,
                                        seg2Node2,
                                        newPos, femSurface,
                                        femSurfaceNorm);
               resetSurfaceProperties = 1;
           }
#endif

           newVel[X] = half * (mergenode1->vX +
                               mergenode2->vX);
           newVel[Y] = half * (mergenode1->vY +
                               mergenode2->vY);
           newVel[Z] = half * (mergenode1->vZ +
                               mergenode2->vZ);


           oldTag1 = mergenode1->myTag;
           oldTag2 = mergenode2->myTag;

           MergeNode(home, OPCLASS_COLLISION, mergenode1,
                     mergenode2, newPos, &targetNode,
                     &mergeStatus, globalOp);

/*
 *         If the merge did not succeed, go back and
 *         continue looking for collision candidates.
 */
           if ((mergeStatus & MERGE_SUCCESS) == 0) {
                mergenode1->flags |= NO_COLLISIONS;
                mergenode2->flags |= NO_COLLISIONS;
                continue;
           }
#ifdef _FEM
/*
 *         Need to explicitly reset surface properties
 *         after colliding 2 surface nodes.
 */
           if (resetSurfaceProperties) {
               targetNode->fem_Surface[0] = femSurface[0];
               targetNode->fem_Surface[1] = femSurface[1];
               targetNode->fem_Surface_Norm[0] =
                       femSurfaceNorm[0];
               targetNode->fem_Surface_Norm[1] =
                       femSurfaceNorm[1];
               targetNode->fem_Surface_Norm[2] =
                       femSurfaceNorm[2];
           }
#endif

/*
 *         When debugging, dump some info on topological
 *         changes taking place and the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
           if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
               printf("  Merge(SegCollision): "
                      "(%d,%d) and (%d,%d) at (%d,%d)\n",
                      oldTag1.domainID, oldTag1.index,
                      oldTag2.domainID, oldTag2.index,
                      targetNode->myTag.domainID,
                      targetNode->myTag.index);
               PrintNode(targetNode);
           }
#endif

/*
 *         If the target node exists after the merge,
 *         reset its velocity, and reset the topological
 *         change exemptions for the target node; use
 *         NodeTopologyExemptions() to get the basic
 *         exemptions, then exempt all the node's arms
 *         from additional segment collisions this cycle.
 */
           if (targetNode != (Node_t *)NULL) {

				if (param->enforceGlidePlanes) {
/*
 *					Handling the collision may have introduced new glide planes
 *					which are being violated. Attempt to fix these violations.
 *					Protect against moving the node too far, however, which
 *					can lead to spurious collisions and nonphysical tangled
 *					dislocation structures.
 */
					real8 oldnewPos[3];
					oldnewPos[X] = newPos[X];
					oldnewPos[Y] = newPos[Y];
					oldnewPos[Z] = newPos[Z];

					AdjustMergePoint(home, targetNode, &newPos[X], &newPos[Y], &newPos[Z]);
					
					vec1[X] = newPos[X] = oldnewPos[X];
					vec1[Y]	= newPos[Y] = oldnewPos[Y];
					vec1[Z]	= newPos[Z] = oldnewPos[Z];
		            if (DotProduct(vec1, vec1) > 16.0 * mindist2) {
        				newPos[X] = oldnewPos[X];
        				newPos[Y] = oldnewPos[Y];
        				newPos[Z] = oldnewPos[Z];
    				}


					FoldBox(param, &newPos[X],&newPos[Y],&newPos[Z]);	

					targetNode->x = newPos[X];
					targetNode->y = newPos[Y];
					targetNode->z = newPos[Z];
				}						

/*
*              If we are enforcing glide planes but
*              allowing some fuzziness in the planes, we
*              also need to recalculate the glide 
*              planes for the segments attched to the
*              collision node.
*/
               if (param->enforceGlidePlanes &&
                   param->allowFuzzyGlidePlanes) {
                   int n;
                   for (n=0; n<targetNode->numNbrs; n++) {
                       tmpNbr = GetNodeFromTag(home,
                               targetNode->nbrTag[n]);
                       RecalcSegGlidePlane(home,
                                           targetNode,
                                           tmpNbr, 1);
                   }
               }
			   
               MarkNodeForceObsolete(home, targetNode);

/*
 *              Set the "older" positions to the new positions for this node
 *              so the collision isn't detected again.
 */
				targetNode->olderx = targetNode->x;
				targetNode->oldery = targetNode->y;
				targetNode->olderz = targetNode->z;

/*
 *             Estimate velocity so mobility function
 *             has a reasonable starting point
 */
               targetNode->vX = newVel[X];
               targetNode->vY = newVel[Y];
               targetNode->vZ = newVel[Z];

               (void)EvaluateMobility(home, targetNode);

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
               targetNode->multiNodeLife = 0;
#endif
           		localCollisionCnt++;

//GenerateOutput(home, STAGE_CYCLE);

           }

/*
 *			Resort the nodes for collisions.
 */
			SortNodes(home,&home->cell2size);

        } 

		} /* while (morecollisions) */

/*
 *		Unflag all nodes - everything should be hunky dory.
 */
    	for (i = 0; i < home->newNodeKeyPtr; i++) {
        	if ((tmpNode = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			tmpNode->flags &= ~NO_COLLISIONS;
		}

/*
 *      Free memory.
 */

        for(i=0;i<MAX_COLS;i++){
			free(done_points[i]);
        }
        free(done_points);
        
        EmptyCollisionList(colList);
		free(colList);
		free(pairlist);

#ifdef DEBUG_LOG_COLLISIONS
#ifdef PARALLEL
        MPI_Reduce(&localCollisionCnt, &globalCollisionCnt, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);
#else
        globalCollisionCnt = localCollisionCnt;
#endif
        if (home->myDomain == 0) {
            printf("  Collision count = %d\n", globalCollisionCnt);
        }
		printf("donesize=%d , colnum=%d\n" , donesize , colnum);
#endif

        TimerStop(home, COLLISION_HANDLING);

        return;
}

#endif
