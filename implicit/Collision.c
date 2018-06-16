/*****************************************************************************
 *
 *      Module:         Collision.c
 *      Description:    This module contains a dispatch function to call the
 *                      user-selected collision handling mechanism plus
 *                      various support functions that are used by all the
 *                      collision handling mechanisms.
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
 *
 *      Included functions:
 *
 *          GetMinDist()
 *          HandleCollisions()
 *
 *****************************************************************************/
#include "Home.h"


#ifdef _FEM
/*---------------------------------------------------------------------------
 *
 *      Function:       FEM_AdjustCollisionPoint
 *      Description:    When two surface nodes collide, the collision point
 *                      must also be on the surface.  This function is used
 *                      to force the collision point to the surface.
 *
 *      Arguments:
 *          seg1Node1,
 *          seg1Node2      Pointers to endpoints of segment 1.  seg1Node1 is
 *                         assumed to be a surface node
 *          seg2Node1    
 *          seg2Node2,     Pointers to endpoints of segment 2.  seg2Node1 is
 *                         assumed to be a surface node
 *          collisionPoint On ebtry to the function contains the collision
 *                         point selected by the caller.  On exit will contain
 *                         a new collision point on the surface.
 *          femSurface     Array in which to return the FEM surface info
 *                         for the surface containing the new collision point
 *          femSurfaceNorm Array in which to return the FEM surface normal
 *                         vector of the surface containing the new collision
 *                         point.
 *
 *-------------------------------------------------------------------------*/
void FEM_AdjustCollisionPoint(Node_t *seg1Node1, Node_t *seg1Node2,
                              Node_t *seg2Node1, Node_t *seg2Node2,
                              real8 collisionPoint[3], int femSurface[2],
                              real8 femSurfaceNorm[3])
{
        real8  lenX, lenY, lenZ, dist1, dist2;
        Node_t *tmpNode;

/*
 *      Assume seg1Node1 and seg2Node1 are both on the surface.  The
 *      new collision point will be the position of the surface node
 *      closest to the originally selected collision point provided
 *      by the caller.
 */
        lenX = seg1Node1->x - collisionPoint[X];
        lenY = seg1Node1->y - collisionPoint[Y];
        lenZ = seg1Node1->z - collisionPoint[Z];

        dist1 = sqrt(lenX * lenX + lenY * lenY + lenZ * lenZ);

        lenX = seg2Node1->x - collisionPoint[X];
        lenY = seg2Node1->y - collisionPoint[Y];
        lenZ = seg2Node1->z - collisionPoint[Z];

        dist2 = sqrt(lenX * lenX + lenY * lenY + lenZ * lenZ);

        tmpNode = ((dist1 < dist2) ? seg1Node1 : seg2Node1);

        collisionPoint[X] = tmpNode->x;
        collisionPoint[Y] = tmpNode->y;
        collisionPoint[Z] = tmpNode->z;

/*
 *      Return the surface properties for the surface containing the
 *      new collision point.
 */
        femSurface[0] = tmpNode->fem_Surface[0];
        femSurface[1] = tmpNode->fem_Surface[1];

        femSurfaceNorm[0] = tmpNode->fem_Surface_Norm[0];
        femSurfaceNorm[1] = tmpNode->fem_Surface_Norm[1];
        femSurfaceNorm[2] = tmpNode->fem_Surface_Norm[2];

        return;
}
#endif


/*---------------------------------------------------------------------------
 *
 *      Function:       GetMinDist
 *      Description:    Find the minimum distance between the two segments
 *                      p1->p2 andp3->p4.
 *
 *      Arguments:
 *          p1x, p1y, p1z  Coordinates of the first endpoint of segment 1
 *          v1x, v1y, v1z  Velocity of the node at point p1
 *          p2x, p2y, p2z  Coordinates of the second endpoint of segment 1
 *          v2x, v2y, v2z  Velocity of the node at point p2
 *          p3x, p3y, p3z  Coordinates of the first endpoint of segment 2
 *          v3x, v3y, v3z  Velocity of the node at point p3
 *          p4x, p4y, p4z  Coordinates of the second endpoint of segment 2
 *          v4x, v4y, v4z  Velocity of the node at point p4
 *          dist2          pointer to location in which to return the square
 *                         of the minimum distance between the two points
 *          ddist2dt       pointter to the location in which to return the
 *                         time rate of change of the distance between
 *                         L1 and L2
 *          L1             pointer to location at which to return the 
 *                         normalized position on seg 1 closest to segment 2
 *          L2             pointer to location at which to return the 
 *                         normalized position on seg 2 closest to segment 1
 *
 *-------------------------------------------------------------------------*/

/***************************************************************************
THIS FUNCTION IS FLAWED AND CAN GIVE INCORRECT SOLUTIONS.  USE GetMinDist2 
BELOW INSTEAD.
****************************************************************************/
void GetMinDist(real8 p1x, real8 p1y, real8 p1z,
                real8 v1x, real8 v1y, real8 v1z,
                real8 p2x, real8 p2y, real8 p2z,
                real8 v2x, real8 v2y, real8 v2z,
                real8 p3x, real8 p3y, real8 p3z,
                real8 v3x, real8 v3y, real8 v3z,
                real8 p4x, real8 p4y, real8 p4z,
                real8 v4x, real8 v4y, real8 v4z,
                real8 *dist2, real8 *ddist2dt, real8 *L1, real8 *L2)
{

      //	Fatal("Using an erroneous version of the GetMinDist function.  Switch the function to GetMinDist2.\n");
        int    i, pos;
        real8  A, B, C, D, E, G;
        real8  eps = 1.0e-12;
        real8  ddistxdt, ddistydt, ddistzdt;
        real8  seg1Lx, seg1Ly, seg1Lz;
        real8  seg2Lx, seg2Ly, seg2Lz;
        real8  seg1vx, seg1vy, seg1vz;
        real8  seg2vx, seg2vy, seg2vz;
        real8  distx, disty, distz;
        real8  dist[4];

        seg1Lx = p2x - p1x;
        seg1Ly = p2y - p1y;
        seg1Lz = p2z - p1z;

        seg2Lx = p4x - p3x;
        seg2Ly = p4y - p3y;
        seg2Lz = p4z - p3z;

        seg1vx = v2x - v1x;
        seg1vy = v2y - v1y;
        seg1vz = v2z - v1z;

        seg2vx = v4x - v3x;
        seg2vy = v4y - v3y;
        seg2vz = v4z - v3z;

        A = seg1Lx*seg1Lx + seg1Ly*seg1Ly + seg1Lz*seg1Lz;
        B = ((seg1Lx*(p1x-p3x)) +
             (seg1Ly*(p1y-p3y)) +
             (seg1Lz*(p1z-p3z))) * 2.0;
        C = (seg1Lx*seg2Lx + seg1Ly*seg2Ly + seg1Lz*seg2Lz) * 2.0;
        D = ((seg2Lx*(p3x-p1x)) +
             (seg2Ly*(p3y-p1y)) +
             (seg2Lz*(p3z-p1z))) * 2.0;
        E = seg2Lx*seg2Lx + seg2Ly*seg2Ly + seg2Lz*seg2Lz;
        G = C*C - 4*A*E;
/*
 *      If segment 1 is just a point...
 */
        if (A < eps) {
            *L1 = 0.0;
            if (E < eps) *L2 = 0.0;
            else *L2 = -0.5 * D / E;

/*
 *      If segment 2 is just a point...
 */
        } else if (E < eps) {
            *L2 = 0.0;
            if (A < eps) *L1 = 0.0;
            else *L1 = -0.5 * B / A;
/*
 *      If segments are parallel
 */
        } else if (fabs(G) < eps) {
            dist[0] = (p3x-p1x)*(p3x-p1x) +
                      (p3y-p1y)*(p3y-p1y) +
                      (p3z-p1z)*(p3z-p1z);
            dist[1] = (p4x-p1x)*(p4x-p1x) +
                      (p4y-p1y)*(p4y-p1y) +
                      (p4z-p1z)*(p4z-p1z);
            dist[2] = (p3x-p2x)*(p3x-p2x) +
                      (p3y-p2y)*(p3y-p2y) +
                      (p3z-p2z)*(p3z-p2z);
            dist[3] = (p4x-p2x)*(p4x-p2x) +
                      (p4y-p2y)*(p4y-p2y) +
                      (p4z-p2z)*(p4z-p2z);

            *dist2 = dist[0];

            for (i = 1; i < 4; i++) {
                if (dist[i] < *dist2) {
                    *dist2 = dist[i];
                    pos = i+1;
                }
            }

            *L1 = floor((real8)pos/2.0);
            *L2 = (real8)(pos-1 % 2);

        } else {
            *L2 = (2.0*A*D + B*C) / G;
            *L1 = 0.5 * (C* *L2 - B) / A;
        }

/*
 *      Make sure L1 and L2 are between 0 and 1
 */
        *L1 = MIN(MAX(*L1, 0.0), 1.0);
        *L2 = MIN(MAX(*L2, 0.0), 1.0);

/*
 *      Lastly, calculate the distance^2 and the time rate of change
 *      between the points at L1 and L2.
 */

	    distx = p1x + (seg1Lx * *L1) - p3x - (seg2Lx * *L2);
	    disty = p1y + (seg1Ly * *L1) - p3y - (seg2Ly * *L2);
	    distz = p1z + (seg1Lz * *L1) - p3z - (seg2Lz * *L2);

	    *dist2 = distx*distx + disty*disty + distz*distz;

        ddistxdt = (v1x + (seg1vx * *L1) - v3x - (seg2vx * *L2)) *
                   (p1x + (seg1Lx * *L1) - p3x - (seg2Lx * *L2));
        ddistydt = (v1y + (seg1vy * *L1) - v3y - (seg2vy * *L2)) *
                   (p1y + (seg1Ly * *L1) - p3y - (seg2Ly * *L2));
        ddistzdt = (v1z + (seg1vz * *L1) - v3z - (seg2vz * *L2)) *
                   (p1z + (seg1Lz * *L1) - p3z - (seg2Lz * *L2));

        *ddist2dt = 2.0 * (ddistxdt + ddistydt + ddistzdt);

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       GetMinDist2
 *      Description:    Find the minimum distance between the two segments
 *                      p1->p2 andp3->p4.
 *
 *      Arguments:
 *          p1x, p1y, p1z  Coordinates of the first endpoint of segment 1
 *          v1x, v1y, v1z  Velocity of the node at point p1
 *          p2x, p2y, p2z  Coordinates of the second endpoint of segment 1
 *          v2x, v2y, v2z  Velocity of the node at point p2
 *          p3x, p3y, p3z  Coordinates of the first endpoint of segment 2
 *          v3x, v3y, v3z  Velocity of the node at point p3
 *          p4x, p4y, p4z  Coordinates of the second endpoint of segment 2
 *          v4x, v4y, v4z  Velocity of the node at point p4
 *          dist2          pointer to location in which to return the square
 *                         of the minimum distance between the two points
 *          ddist2dt       pointter to the location in which to return the
 *                         time rate of change of the distance between
 *                         L1 and L2
 *          L1             pointer to location at which to return the 
 *                         normalized position on seg 1 closest to segment 2
 *          L2             pointer to location at which to return the 
 *                         normalized position on seg 2 closest to segment 1
 *
 *-------------------------------------------------------------------------*/

void GetMinDist2(real8 p1x, real8 p1y, real8 p1z,
                real8 v1x, real8 v1y, real8 v1z,
                real8 p2x, real8 p2y, real8 p2z,
                real8 v2x, real8 v2y, real8 v2z,
                real8 p3x, real8 p3y, real8 p3z,
                real8 v3x, real8 v3y, real8 v3z,
                real8 p4x, real8 p4y, real8 p4z,
                real8 v4x, real8 v4y, real8 v4z,
                real8 *dist2, real8 *ddist2dt, real8 *L1, real8 *L2)
{
        int    i, pos;
		int 	icase, didDist2;
        real8  A, B, C, D, E, G;
        real8  eps = 1.0e-12;
        real8  ddistxdt, ddistydt, ddistzdt;
        real8  seg1Lx, seg1Ly, seg1Lz;
        real8  seg2Lx, seg2Ly, seg2Lz;
        real8  seg1vx, seg1vy, seg1vz;
        real8  seg2vx, seg2vy, seg2vz;
        real8  distx, disty, distz, d2, d2min;
        real8  dist[4];

		real8	r1mr3[3], r2mr1[3], r4mr3[3], r4mr1[3], r3mr2[3], r4mr2[3];
		real8	M[2][2], rhs[2], sol[2], detM;
		real8   trial[4][2];

		r1mr3[0] = p1x-p3x; r1mr3[1] = p1y-p3y; r1mr3[2] = p1z-p3z; 
		r2mr1[0] = p2x-p1x; r2mr1[1] = p2y-p1y; r2mr1[2] = p2z-p1z; 
		r4mr3[0] = p4x-p3x; r4mr3[1] = p4y-p3y; r4mr3[2] = p4z-p3z;
		
		seg1Lx = r2mr1[0] ;  seg1Ly = r2mr1[1] ;  seg1Lz = r2mr1[2];
        seg2Lx = r4mr3[0] ;  seg2Ly = r4mr3[1] ;  seg2Lz = r4mr3[2];

        seg1vx = v2x - v1x;  seg1vy = v2y - v1y; seg1vz = v2z - v1z;
        seg2vx = v4x - v3x;  seg2vy = v4y - v3y; seg2vz = v4z - v3z;

		M[0][0] = DotProduct(r2mr1,r2mr1);
		M[1][0] =-DotProduct(r4mr3,r2mr1);
		M[1][1] = DotProduct(r4mr3,r4mr3);
		M[0][1] = M[1][0];

		rhs[0] = -DotProduct(r2mr1,r1mr3);
		rhs[1] =  DotProduct(r4mr3,r1mr3);

	//	detM = M[0][0]*M[1][1] - M[0][1]*M[1][0];
		detM = 1.0 - M[1][0] * M[1][0] / M[0][0] / M[1][1];
		
	//	if (detM<0) printf("Warning: negative detM = %e in GetMinDist2\n" , detM);
	//	if (detM<0) detM=eps/2.0;

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
            *L1 = 0.0;
            if (E < eps) *L2 = 0.0;
            else *L2 = -0.5 * D / E;

/*
 *      If segment 2 is just a point...
 */
        } else if (E < eps) {
            *L2 = 0.0;
            if (A < eps) *L1 = 0.0;
            else *L1 = -0.5 * B / A;
/*
 *      If segments are parallel
 */
		} else if (detM<1e-6) {
            r4mr1[0] = p4x-p1x;  r4mr1[1] = p4y-p1y;  r4mr1[2] = p4z-p1z;
			r3mr2[0] = p3x-p2x;  r3mr2[1] = p3y-p2y;  r3mr2[2] = p3z-p2z;
			r4mr2[0] = p4x-p2x;  r4mr2[1] = p4y-p2y;  r4mr2[2] = p4z-p2z;
			
			dist[0] = DotProduct(r1mr3,r1mr3);
            dist[1] = DotProduct(r4mr1,r4mr1);
            dist[2] = DotProduct(r3mr2,r3mr2);
            dist[3] = DotProduct(r4mr2,r4mr2);

            *dist2 = dist[0];
            pos = 1;

            for (i = 1; i < 4; i++) {
                if (dist[i] < *dist2) {
                    *dist2 = dist[i];
                    pos = i+1;
                }
            }

            *L1 = floor((real8)pos/2.1);
			*L2 = (real8)(1 - (pos % 2));
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
				*L1 = sol[0];
				*L2 = sol[1];

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
//				printf("%e %e\n%e %e\n%e %e\n%e %e\n",trial[0][0],trial[0][1],
	//						trial[1][0],trial[1][1],trial[2][0],trial[2][1],
		//					trial[3][0],trial[3][1]);

				/* find the minimum out of four trials */
				d2min = 1e100;
				for(icase = 0; icase < 4; icase++) {
					trial[icase][0] = MIN(MAX(trial[icase][0], 0.0), 1.0);
					trial[icase][1] = MIN(MAX(trial[icase][1], 0.0), 1.0);
					distx = p1x + (seg1Lx * trial[icase][0]) 
					      - p3x - (seg2Lx * trial[icase][1]);
					disty = p1y + (seg1Ly * trial[icase][0]) 
					      - p3y - (seg2Ly * trial[icase][1]);
					distz = p1z + (seg1Lz * trial[icase][0]) 
					      - p3z - (seg2Lz * trial[icase][1]);

					d2 = distx*distx + disty*disty + distz*distz;
					if (d2<d2min) {
						*L1 = trial[icase][0];
						*L2 = trial[icase][1];
						d2min = d2;
					}
				}
				*dist2 = d2min;
				didDist2 = 1;
			}
		} 

/*
 *      Make sure L1 and L2 are between 0 and 1
 */
        *L1 = MIN(MAX(*L1, 0.0), 1.0);
        *L2 = MIN(MAX(*L2, 0.0), 1.0);

		if (!didDist2) {
			distx = p1x + (seg1Lx * *L1) - p3x - (seg2Lx * *L2);
			disty = p1y + (seg1Ly * *L1) - p3y - (seg2Ly * *L2);
			distz = p1z + (seg1Lz * *L1) - p3z - (seg2Lz * *L2);

			*dist2 = distx*distx + disty*disty + distz*distz;
		}

        ddistxdt = (v1x + (seg1vx * *L1) - v3x - (seg2vx * *L2)) *
                   (p1x + (seg1Lx * *L1) - p3x - (seg2Lx * *L2));
        ddistydt = (v1y + (seg1vy * *L1) - v3y - (seg2vy * *L2)) *
                   (p1y + (seg1Ly * *L1) - p3y - (seg2Ly * *L2));
        ddistzdt = (v1z + (seg1vz * *L1) - v3z - (seg2vz * *L2)) *
                   (p1z + (seg1Lz * *L1) - p3z - (seg2Lz * *L2));

        *ddist2dt = 2.0 * (ddistxdt + ddistydt + ddistzdt);

        return;
}


void HandleCollisions(Home_t *home)
{
	int passInt;
        switch(home->param->collisionMethod) {
            case 3:
                RetroactiveCollisions(home);
                break;

            case 2:
                PredictiveCollisions(home);
                break;

            case 1:
            default:
                ProximityCollisions(home);
                break;
        }

        return;
}
