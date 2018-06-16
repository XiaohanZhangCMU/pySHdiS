#include "Home.h"
#include "Util.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#define SEG_NUM_VALS 13

int isCollinear(double nX, double nY, double nZ, double pX, double pY, double pZ);
int GetSlipSystem(Home_t *home, double bX, double bY, double bZ, double nX, double nY, double nZ);
int OutsideVolume(Home_t *home, double px, double py, double pz);
int FacetIntersectionPosition(Home_t *home, double xin, double yin, double zin, 
                              double xout, double yout, double zout, 
					          double *xpos, double *ypos, double *zpos);

/*---------------------------------------------------------------------------
 *
 *      Function:    WriteParaview
 *      Description: Write segment data in ParaView (vtk) format (Nicolas)
 *
 *      Args:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *
 *-------------------------------------------------------------------------*/
void WriteParaview(Home_t *home, char *baseFileName, int numSegs, int numLocSegs)
{
        int      i, j, thisDomain;
        int      newNodeKeyPtr, narm;
        int      btype, out, facet;
        real8    x1, y1, z1;
        real8    x2, y2, z2;
        real8    d1x, d1y, d1z;
        real8    d2x, d2y, d2z;
        real8    vx, vy, vz, vmag;
        real8    xpos, ypos, zpos;
        real8    bX, bY, bZ;
        real8    nX, nY, nZ;
        real8    Lx, Ly, Lz;
        char     fileName[256];
        Node_t   *node, *nbrNode;
        Param_t  *param;
        FILE     *fp;
        struct   stat statbuf;
        int      nLocalSeg, maxLocSegs, nGlobalSegs;
        double   *localListSeg;
        double   *globalListSeg;
#ifdef PARALLEL
        MPI_Status status;
#endif


        param      = home->param;
        thisDomain = home->myDomain;

        Lx = param->Lx;
        Ly = param->Ly;
        Lz = param->Lz;

/*
 *      Set data file name.
 */
        snprintf(fileName, sizeof(fileName), "%s/%s.vtk",
                 DIR_PARAVIEW, baseFileName);


/*
 *      Only the first task (myDomain=0) outputs ParaView files
 */
        if (thisDomain == 0) {
/*
 *          First task must open the data file for writing
 *          to overwrite any existing file of the same name.
 */
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("Paraview: Open error %d on %s\n", errno, fileName);
            }

			printf(" +++ Writing Paraview file(s) %s\n", baseFileName);
			
			fprintf(fp, "# vtk DataFile Version 3.0\n");
			fprintf(fp, "ParaDiS dislocation segments\n");
			fprintf(fp, "ASCII\n");
			fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

        }
 

/*
 *      Generate the plot data for the segments associated with this
 *      domain's data
 */
        maxLocSegs = 2*numLocSegs;
        localListSeg = (double*)malloc(sizeof(double)*SEG_NUM_VALS*maxLocSegs);
        nLocalSeg = 0;
        
        newNodeKeyPtr = home->newNodeKeyPtr;

        for (i = 0; i < newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {;
                continue;
            }
        
            x1 = node->x;
            y1 = node->y;
            z1 = node->z;
            
            d1x = x1 - node->oldx;
			d1y = y1 - node->oldy;
			d1z = z1 - node->oldz;
			ZImage(home->param, &d1x, &d1y, &d1z);
        
            for (j = 0; j < node->numNbrs; j++) {
        
/*
 *  FIX ME? This will result in segments crossing domain boundaries
 *  to be added to the file by both domains!  Is this what is wanted?
 */
                if ((node->nbrTag[j].domainID == thisDomain) && 
                    (node->nbrTag[j].index < i)) {
                    continue;
                }

                bX = node->burgX[j];
                bY = node->burgY[j];
                bZ = node->burgZ[j];
			
                nX = node->nx[j];
                nY = node->ny[j];
                nZ = node->nz[j];
                
/*
 *              For the following burgers vector checks, convert the burgers
 *              vector to the crystalographic frame if necessary.
 */
                if (param->useLabFrame) {
                    real8 burgLab[3] = {bX, bY, bZ}, burgCrystal[3];
                    real8 normLab[3] = {nX, nY, nZ}, normCrystal[3];

                    Matrix33Vector3Multiply(home->rotMatrixInverse, burgLab,
                                            burgCrystal);
                    bX = burgCrystal[X];
                    bY = burgCrystal[Y];
                    bZ = burgCrystal[Z];
                    
                    Matrix33Vector3Multiply(home->rotMatrixInverse, normLab,
                                            normCrystal);
                    nX = normCrystal[X];
                    nY = normCrystal[Y];
                    nZ = normCrystal[Z];
                }
        
                nbrNode = GetNeighborNode(home, node, j);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }
                
                x2 = nbrNode->x;
                y2 = nbrNode->y;
                z2 = nbrNode->z;
                
                d2x = x2 - nbrNode->oldx;
				d2y = y2 - nbrNode->oldy;
				d2z = z2 - nbrNode->oldz;
				ZImage(home->param, &d2x, &d2y, &d2z);
				
				vx = 0.5*(d1x + d2x)/param->realdt;
				vy = 0.5*(d1y + d2y)/param->realdt;
				vz = 0.5*(d1z + d2z)/param->realdt;
				vmag = sqrt(vx*vx + vy*vy + vz*vz);
                
                PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
                
                out = OutsideVolume(home, x2, y2, z2);
                while (out) {
					facet = FacetIntersectionPosition(home, x1, y1, z1, x2, y2, z2, 
					                                  &xpos, &ypos, &zpos);
					 
					if (facet == -1) {
						out = 0;
						continue;
					}
					
					// Add segment
					if (nLocalSeg >= maxLocSegs) {
						maxLocSegs += 100;
						localListSeg = (double*)realloc(localListSeg, 
						               sizeof(double)*SEG_NUM_VALS*maxLocSegs);
						
						/* Just to make sure it doesn't loop indefinitely... */
						if (maxLocSegs >= 100*numLocSegs) {
							Fatal("Unexpected number of segments in WriteParaview");
						}
					}
					
					localListSeg[nLocalSeg*SEG_NUM_VALS+0] = x1;
					localListSeg[nLocalSeg*SEG_NUM_VALS+1] = y1;
					localListSeg[nLocalSeg*SEG_NUM_VALS+2] = z1;
					localListSeg[nLocalSeg*SEG_NUM_VALS+3] = xpos;
					localListSeg[nLocalSeg*SEG_NUM_VALS+4] = ypos;
					localListSeg[nLocalSeg*SEG_NUM_VALS+5] = zpos;
					localListSeg[nLocalSeg*SEG_NUM_VALS+6] = bX;
					localListSeg[nLocalSeg*SEG_NUM_VALS+7] = bY;
					localListSeg[nLocalSeg*SEG_NUM_VALS+8] = bZ;
					localListSeg[nLocalSeg*SEG_NUM_VALS+9] = nX;
					localListSeg[nLocalSeg*SEG_NUM_VALS+10] = nY;
					localListSeg[nLocalSeg*SEG_NUM_VALS+11] = nZ;
					localListSeg[nLocalSeg*SEG_NUM_VALS+12] = vmag;
					nLocalSeg++;
					
					if      (facet == 0) xpos = param->maxSideX;
					else if (facet == 1) ypos = param->maxSideY;
					else if (facet == 2) zpos = param->maxSideZ;
					else if (facet == 3) xpos = param->minSideX;
					else if (facet == 4) ypos = param->minSideY;
					else if (facet == 5) zpos = param->minSideZ;
					
					x1 = xpos;
					y1 = ypos;
					z1 = zpos;
					
					x2 = nbrNode->x;
					y2 = nbrNode->y;
					z2 = nbrNode->z;
					
					PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
					out = OutsideVolume(home, x2, y2, z2);
				}
/*
 *              Add this segment to the local domain's list 
 */
                if (nLocalSeg >= maxLocSegs) {
					maxLocSegs += 100;
					localListSeg = (double*)realloc(localListSeg, 
					               sizeof(double)*SEG_NUM_VALS*maxLocSegs);
				}
                
                localListSeg[nLocalSeg*SEG_NUM_VALS+0] = x1;
                localListSeg[nLocalSeg*SEG_NUM_VALS+1] = y1;
                localListSeg[nLocalSeg*SEG_NUM_VALS+2] = z1;
                localListSeg[nLocalSeg*SEG_NUM_VALS+3] = x2;
                localListSeg[nLocalSeg*SEG_NUM_VALS+4] = y2;
                localListSeg[nLocalSeg*SEG_NUM_VALS+5] = z2;
                localListSeg[nLocalSeg*SEG_NUM_VALS+6] = bX;
                localListSeg[nLocalSeg*SEG_NUM_VALS+7] = bY;
                localListSeg[nLocalSeg*SEG_NUM_VALS+8] = bZ;
                localListSeg[nLocalSeg*SEG_NUM_VALS+9] = nX;
                localListSeg[nLocalSeg*SEG_NUM_VALS+10] = nY;
                localListSeg[nLocalSeg*SEG_NUM_VALS+11] = nZ;
                localListSeg[nLocalSeg*SEG_NUM_VALS+12] = vmag;
                nLocalSeg++;
                
            }
        }
        

#ifdef PARALLEL
		MPI_Reduce(&nLocalSeg, &nGlobalSegs, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);
#else
		nGlobalSegs = nLocalSeg;
#endif       
        if (thisDomain == 0) {
			globalListSeg = (double*)malloc(sizeof(double)*SEG_NUM_VALS*nGlobalSegs);
		}
        
#ifdef PARALLEL
		if (thisDomain == 0) {
			
			for(i = 0; i < SEG_NUM_VALS*nLocalSeg; i++) {
				globalListSeg[i] = localListSeg[i];
			}
			int indSeg = SEG_NUM_VALS*nLocalSeg;
			
			// Receive distant segments from other tasks
			for(i = 1; i < home->numDomains; i++) {
				int numDistSegs = -1;
				MPI_Recv(&numDistSegs, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
				if (numDistSegs > 0) {
					double *distListSeg = (double*)malloc(sizeof(double)*SEG_NUM_VALS*numDistSegs);
					MPI_Recv(distListSeg, SEG_NUM_VALS*numDistSegs, MPI_DOUBLE, 
					         status.MPI_SOURCE, 2, MPI_COMM_WORLD, &status);
					for(j = 0; j < SEG_NUM_VALS*numDistSegs; j++) {
						globalListSeg[indSeg++] = distListSeg[j];
					}
					free(distListSeg);
				}
			}
		}
		else {
			// Send local segment to first task
			MPI_Send(&nLocalSeg, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
			if (nLocalSeg > 0) {
				MPI_Send(localListSeg, SEG_NUM_VALS*nLocalSeg, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
			}
		}
#else
		for(i = 0; i < SEG_NUM_VALS*nLocalSeg; i++) {
			globalListSeg[i] = localListSeg[i];
		}
#endif
        
        free(localListSeg);

/*
 *      Only the first task (myDomain=0) outputs ParaView files
 */
        if (thisDomain == 0) {
			
			fprintf(fp, "POINTS %d float\n", 2*nGlobalSegs + 8);
	
			// Simulation box vertices
			fprintf(fp, "%f %f %f\n", param->minSideX, param->minSideY, param->minSideZ);
			fprintf(fp, "%f %f %f\n", param->maxSideX, param->minSideY, param->minSideZ);
			fprintf(fp, "%f %f %f\n", param->maxSideX, param->maxSideY, param->minSideZ);
			fprintf(fp, "%f %f %f\n", param->minSideX, param->maxSideY, param->minSideZ);
			fprintf(fp, "%f %f %f\n", param->minSideX, param->minSideY, param->maxSideZ);
			fprintf(fp, "%f %f %f\n", param->maxSideX, param->minSideY, param->maxSideZ);
			fprintf(fp, "%f %f %f\n", param->maxSideX, param->maxSideY, param->maxSideZ);
			fprintf(fp, "%f %f %f\n", param->minSideX, param->maxSideY, param->maxSideZ);
	
			// Dislocation nodes position
			for(i = 0; i < nGlobalSegs; i++) {
				fprintf(fp, "%f %f %f\n", globalListSeg[i*SEG_NUM_VALS+0],
				globalListSeg[i*SEG_NUM_VALS+1], globalListSeg[i*SEG_NUM_VALS+2]);
				fprintf(fp, "%f %f %f\n", globalListSeg[i*SEG_NUM_VALS+3],
				globalListSeg[i*SEG_NUM_VALS+4], globalListSeg[i*SEG_NUM_VALS+5]);
			}
	
			// Drawing the simulation box
			fprintf(fp, "CELLS %d %d\n", 1 + nGlobalSegs, 9 + 3 * nGlobalSegs);
			fprintf(fp, "8 0 1 2 3 4 5 6 7\n");
	
			// Drawing the dislocation segments
			for(i = 0; i < nGlobalSegs; i++) {
				fprintf(fp, "%d %d %d\n", 2, 2*i + 8, 2*i + 9);
			}

			// Cells type for simulation box
			fprintf(fp, "CELL_TYPES %d\n", 1 + nGlobalSegs);
			fprintf(fp, "12\n");
	
			// Cells type for dislocation segments
			for(i = 0; i < nGlobalSegs; i++) {
				fprintf(fp, "4\n");
			}
			
			// Cell data
			fprintf(fp, "CELL_DATA %d\n", nGlobalSegs + 1);
			fprintf(fp, "SCALARS Slip_System int 1\n");
			fprintf(fp, "LOOKUP_TABLE default\n");
			fprintf(fp, "-1\n");
			for(i = 0; i < nGlobalSegs; i++) {
				int sys;
				bX = globalListSeg[i*SEG_NUM_VALS+6];
				bY = globalListSeg[i*SEG_NUM_VALS+7];
				bZ = globalListSeg[i*SEG_NUM_VALS+8];
				nX = globalListSeg[i*SEG_NUM_VALS+9];
				nY = globalListSeg[i*SEG_NUM_VALS+10];
				nZ = globalListSeg[i*SEG_NUM_VALS+11];
				sys = GetSlipSystem(home, bX, bY, bZ, nX, nY, nZ);
				fprintf(fp, "%d\n", sys);
			}
			fprintf(fp, "SCALARS Junction_type int 1\n");
			fprintf(fp, "LOOKUP_TABLE default\n");
			fprintf(fp, "-1\n");
			for(i = 0; i < nGlobalSegs; i++) {
				int type = -1;
				nX = globalListSeg[i*SEG_NUM_VALS+9];
				nY = globalListSeg[i*SEG_NUM_VALS+10];
				nZ = globalListSeg[i*SEG_NUM_VALS+11];
				if (home->param->materialType == MAT_TYPE_FCC) {
					if (isCollinear(nX, nY, nZ, 1, 1, 1) ||
					    isCollinear(nX, nY, nZ, -1, 1, 1) ||
					    isCollinear(nX, nY, nZ, 1, -1, 1) ||
					    isCollinear(nX, nY, nZ, 1, 1, -1)) {
						type = 0; // Glissile segments
					}
					else {
						bX = globalListSeg[i*SEG_NUM_VALS+6];
						bY = globalListSeg[i*SEG_NUM_VALS+7];
						bZ = globalListSeg[i*SEG_NUM_VALS+8];
						if (isCollinear(nX, nY, nZ, 1, 0, 0) ||
						    isCollinear(nX, nY, nZ, 0, 1, 0) ||
						    isCollinear(nX, nY, nZ, 0, 0, 1)) {
							type = 1; // Lomer
						}
						else if (isCollinear(bX, bY, bZ, 1, 0, 0) ||
						         isCollinear(bX, bY, bZ, 0, 1, 0) ||
						         isCollinear(bX, bY, bZ, 0, 0, 1)) {
							type = 2; // Hirth
						}
						else {
							type = 3; // Unidentified junctions
						}
					}
				}
				else {
					type = 0;
				}
				fprintf(fp, "%d\n", type);
			}
			fprintf(fp, "VECTORS Burgers FLOAT\n");
			fprintf(fp, "0.0 0.0 0.0\n");
			for(i = 0; i < nGlobalSegs; i++) {
				fprintf(fp, "%f %f %f\n", globalListSeg[i*SEG_NUM_VALS+6], 
				globalListSeg[i*SEG_NUM_VALS+7], globalListSeg[i*SEG_NUM_VALS+8]);
			}
			fprintf(fp, "NORMALS Normals FLOAT\n");
			fprintf(fp, "0.0 0.0 0.0\n");
			for(i = 0; i < nGlobalSegs; i++) {
				fprintf(fp, "%f %f %f\n", globalListSeg[i*SEG_NUM_VALS+9], 
				globalListSeg[i*SEG_NUM_VALS+10], globalListSeg[i*SEG_NUM_VALS+11]);
			}
			fprintf(fp, "SCALARS Velocity FLOAT\n");
			fprintf(fp, "LOOKUP_TABLE default\n");
			fprintf(fp, "0.0\n");
			for(i = 0; i < nGlobalSegs; i++) {
				fprintf(fp, "%f\n", globalListSeg[i*SEG_NUM_VALS+12]);
			}
			
            fclose(fp);
            free(globalListSeg);
        }
       
        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:    isCollinear
 *
 *-------------------------------------------------------------------------*/
int isCollinear(double nX, double nY, double nZ, double pX, double pY, double pZ)
{
	double n1[3], n2[3], n1n2;
	
	n1[0] = nX;
	n1[1] = nY;
	n1[2] = nZ;
	NormalizeVec(n1);
	
	n2[0] = pX;
	n2[1] = pY;
	n2[2] = pZ;
	NormalizeVec(n2);
	
	n1n2 = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
	if (fabs(fabs(n1n2)-1.0) < 1.0e-5) return 1;
	
	return 0;
}

/*---------------------------------------------------------------------------
 *
 *      Function:    GetSlipSystem
 *
 *-------------------------------------------------------------------------*/
int GetSlipSystem(Home_t *home, double bX, double bY, double bZ, 
                                double nX, double nY, double nZ)
{
	int nid, bid, sys;
	
	if (home->param->materialType == MAT_TYPE_FCC) {
		nid = -1;
		if      (isCollinear(nX, nY, nZ, -1.0, 1.0, 1.0))  nid = 0;
		else if (isCollinear(nX, nY, nZ, 1.0, 1.0, 1.0))   nid = 1;
		else if (isCollinear(nX, nY, nZ, -1.0, -1.0, 1.0)) nid = 2;
		else if (isCollinear(nX, nY, nZ, 1.0, -1.0, 1.0))  nid = 3;
		bid = -1;
		if      (isCollinear(bX, bY, bZ, 0.0, 1.0, 1.0))  bid = 0;
		else if (isCollinear(bX, bY, bZ, 0.0, -1.0, 1.0)) bid = 1;
		else if (isCollinear(bX, bY, bZ, 1.0, 0.0, 1.0))  bid = 2;
		else if (isCollinear(bX, bY, bZ, -1.0, 0.0, 1.0)) bid = 3;
		else if (isCollinear(bX, bY, bZ, -1.0, 1.0, 0.0)) bid = 4;
		else if (isCollinear(bX, bY, bZ, 1.0, 1.0, 0.0))  bid = 5;

		if (nid == 0) {
			if      (bid == 1) sys = 1; // A2
			else if (bid == 2) sys = 2; // A3
			else if (bid == 5) sys = 3; // A6
			else sys = 0;
		} else if (nid == 1) {
			if      (bid == 1) sys = 4; // B2
			else if (bid == 3) sys = 5; // B4
			else if (bid == 4) sys = 6; // B5
			else sys = 0;
		} else if (nid == 2) {
			if      (bid == 0) sys = 7; // C1
			else if (bid == 2) sys = 8; // C3
			else if (bid == 4) sys = 9; // C5
			else sys = 0;
		} else if (nid == 3) {
			if      (bid == 0) sys = 10; // D1
			else if (bid == 3) sys = 11; // D4
			else if (bid == 5) sys = 12; // D6
			else sys = 0;
		} else {
			sys = 0; // Junction or non-native FCC system
		}
	} else {
		sys = 0;
	}
	
	return sys;
}

/*---------------------------------------------------------------------------
 *
 *      Function:    OutsideVolume
 *
 *-------------------------------------------------------------------------*/
int OutsideVolume(Home_t *home, double px, double py, double pz)
{
	Param_t *param;
	int out = 0;
	
	param = home->param;
	
	if        (param->xBoundType == Periodic && 
	           (px < param->minSideX || px > param->maxSideX)) {
		out = 1;
	} else if (param->yBoundType == Periodic && 
	           (py < param->minSideY || py > param->maxSideY)) {
		out = 1;
	} else if (param->zBoundType == Periodic && 
	           (pz < param->minSideZ || pz > param->maxSideZ)) {
		out = 1;
	}
	
	return out;
}

/*---------------------------------------------------------------------------
 *
 *      Function:    FacetIntersectionPosition
 *
 *-------------------------------------------------------------------------*/
int FacetIntersectionPosition(Home_t *home, double xin, double yin, double zin, 
                              double xout, double yout, double zout, 
					          double *xpos, double *ypos, double *zpos)
{
	int     i;
	real8   bounds[6];
	real8   p[3], n[3], w[3], t[3], D, N, st;
	int     fint[6];
	Param_t *param;
	
	int facet = -1;
	
	*xpos = xout;
	*ypos = yout;
	*zpos = zout;
	
	t[0] = xout - xin;
	t[1] = yout - yin;
	t[2] = zout - zin;
	
	param = home->param;
	
	bounds[0] = param->minSideX;
	bounds[1] = param->minSideY;
	bounds[2] = param->minSideZ;
	bounds[3] = param->maxSideX;
	bounds[4] = param->maxSideY;
	bounds[5] = param->maxSideZ;
	
	fint[0] = (xout < bounds[0]);
	fint[1] = (yout < bounds[1]);
	fint[2] = (zout < bounds[2]);
	fint[3] = (xout > bounds[3]);
	fint[4] = (yout > bounds[4]);
	fint[5] = (zout > bounds[5]);
	
	for (i = 0; i < 6; i++) {
		if (!fint[i]) continue;
		
		// Point on the surface
		p[0] = 0.0;
		p[1] = 0.0;
		p[2] = 0.0;
		
		// Normal to the surface
		n[0] = 0.0;
		n[1] = 0.0;
		n[2] = 0.0;
		
		if (i < 3) {
			p[i] = bounds[i];
			n[i] = -1.0;
		} else {
			p[i-3] = bounds[i];
			n[i-3] = 1.0;
		}
		
		w[0] = xin - p[0];
		w[1] = yin - p[1];
		w[2] = zin - p[2];
		
		D =  n[0]*t[0] + n[1]*t[1] + n[2]*t[2];
		N = -n[0]*w[0] - n[1]*w[1] - n[2]*w[2];
		
		// Segment is parallel to the plane
		if (fabs(D) < 1.e-5) {
			continue; // No intersection
		}
		
		// Compute intersection parameter
		st = N / D;
		if (st < 0.0 || st > 1.0) continue;
		
		// Compute intersection point
		*xpos = xin + st * t[0];
		*ypos = yin + st * t[1];
		*zpos = zin + st * t[2];
		
		if (!OutsideVolume(home, *xpos, *ypos, *zpos)) {
			facet = i;
			break;
		} else {
			*xpos = xout;
			*ypos = yout;
			*zpos = zout;
		}
	}
	
	return facet;
}
