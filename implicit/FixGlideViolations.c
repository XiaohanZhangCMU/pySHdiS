/*------------------------------------------------------------------------
 *
 *      Function:    FixGlideViolations
 *      Description: Orthogonalize current nodal positions onto glide planes
 *		defined by old nodal positions. Designed to correct
 *		glide plane violations as they occur and prevent them from
 *		accumulating.
 *
 *      Two modes of operation:
 *      1) Provide a node tag, and only correct that node. If an old position
 *         is passed in also oldpos, use that but otherwise use the node's
 *         old position.
 *      2) A NULL tag is provided, try to fix violations for all nodes.
 *
 *		Witten by: Ryan Sills, 4/23/14
 *
 *-----------------------------------------------------------------------*/
#include "Home.h"


void FixGlideViolations(Home_t *home, Tag_t *tag, real8 oldpos[3]) 
{

	int		i, k, l, numplanes, imax;
	real8	nxtmp, nytmp, nztmp;
	real8	x, y, z, amin, eqn, d, nx, ny, nz;
    real8   oldx, oldy, oldz;
	real8	planes[100]; //arbitrarily large
	real8	line[3], dr[3], tmpvec1[3], tmpvec2[3];
	real8 	tol = 1e-5;
	Node_t	*node, *nbr;
	Param_t	*param;

	param = home->param;

	if (param->enforceGlidePlanes) {
        if (tag == NULL) {
            imax = home->newNodeKeyPtr;
        } else {
            imax = 1;
        }
		for (i = 0; i < imax; i++) {

            if (tag == NULL) {
			    if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
                oldx = node->oldx;
                oldy = node->oldy;
                oldz = node->oldz;
            } else {
                node = GetNodeFromTag (home, *tag);
                if (oldpos == NULL) {
                    oldx = node->oldx;
                    oldy = node->oldy;
                    oldz = node->oldz;                   
                } else {
                    oldx = oldpos[X];
                    oldy = oldpos[Y];
                    oldz = oldpos[Z];
                }
            }
/*
 *			If this nodes is fixed we can't move it.
 */
			if (node->constraint == PINNED_NODE) continue;

/*
 *			Determine the number of unique glide planes.
 */
			numplanes = 0;
			for (k = 0; k < node->numNbrs; k++) {
				nxtmp = node->nx[k];
				nytmp = node->ny[k];
				nztmp = node->nz[k];
				Normalize(&nxtmp,&nytmp,&nztmp);
				if (k==0) {
					planes[0] = nxtmp;
					planes[1] = nytmp;
					planes[2] = nztmp;
					numplanes = 1;
				} else {
					for (l=0; l<numplanes; l++) {
						Orthogonalize(&nxtmp,&nytmp,&nztmp,
										planes[3*l],planes[3*l+1],planes[3*l+2]);
					}
					if ((nxtmp*nxtmp+nytmp*nytmp+nztmp*nztmp)>tol) {
						Normalize(&nxtmp,&nytmp,&nztmp);
						planes[3*numplanes] = nxtmp;
						planes[3*numplanes+1] = nytmp;
						planes[3*numplanes+2] = nztmp;
						numplanes++;
					}
				}
			}

/*
 *			Now we orthogonalize. If the node is only on one plane, we project
 *			it onto that plane. If it is on two planes, we project it onto
 *			the line defined by their intersection. If it is on more than two
 *			planes, it should not have moved so reposition it back to its old
 *			position.
 */
			x = node->x; y = node->y; z = node->z;			
			PBCPOSITION(param,oldx,oldy,oldz,&x,&y,&z);

			if (numplanes==1) {
				nx = planes[0]; ny = planes[1]; nz = planes[2];

				d = nx*oldx + ny*oldy + nz*oldz;

				eqn = nx*x + ny*y + nz*z - d; 

				x = x - nx*eqn;
				y = y - ny*eqn;
				z = z - nz*eqn;
			} else if (numplanes==2) {
				tmpvec1[0] = planes[0]; tmpvec1[1] = planes[1]; tmpvec1[2] = planes[2];
				tmpvec2[0] = planes[3]; tmpvec2[1] = planes[4]; tmpvec2[2] = planes[5];
				NormalizedCrossVector(tmpvec1, tmpvec2, line);	

				dr[0] = x-oldx;
				dr[1] = y-oldy;	
				dr[2] = z-oldz;

				amin = DotProduct(line,dr);

				x = oldx + amin*line[0];
				y = oldy + amin*line[1];
				z = oldz + amin*line[2];
			} else {
				x = oldx;
				y = oldy;
				z = oldz;
			}

			FoldBox(param,&x,&y,&z);

			node->x = x;
			node->y = y;
			node->z = z;
		}
	}
	return;
}

/*---------------------------------------------------------------------------
 *
 *	Function:	    TestForGlideViolations
 *	Description:	Loop over all nodes and look for glide plane violations.
 *                  Used for debugging.
 *
 *-------------------------------------------------------------------------*/

void TestForGlideViolations(Home_t *home, char loc_string[]) {

		int		i, j;
		real8	nx, ny, nz, nbrx, nbry, nbrz, drx, dry, drz, drlen;
		real8	tol_test, tol_len, test_len, dottest, devlen;
		Node_t	*node, *nbr;
		Param_t	*param;

		param = home->param;
		tol_test = 1e-5;  //tolerance for plane test, fraction of rann
		tol_len = 1e-10;  //tolerance for "zero-length" segment

        if (param->enforceGlidePlanes) {

            test_len = tol_test*param->rann;
		    for (i = 0; i < home->newNodeKeyPtr; i++) {
			    if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
		        for (j=0; j < node->numNbrs; j++) {

			        nx = node->nx[j];
			        ny = node->ny[j];
			        nz = node->nz[j];
			        Normalize(&nx,&ny,&nz);

			        nbr = GetNeighborNode(home, node, j);

					//if (OrderNodes(node, nbr) >= 0) continue;

			        drx = node->x - nbr->x;
			        dry = node->y - nbr->y;
			        drz = node->z - nbr->z;
                    ZImage(param,&drx,&dry,&drz);

			        drlen = sqrt(drx*drx+dry*dry+drz*drz);
			        if (drlen<tol_len) continue;

			        Normalize(&drx,&dry,&drz);
			        dottest = fabs(nx*drx+ny*dry+nz*drz);

                    devlen = dottest*drlen; 

			        if (devlen > test_len) {
				        printf("WARNING: Node %i,%i glide plane violation after %s " 
                               "out of plane by height %e.\n",
                               node->myTag.domainID,node->myTag.index,loc_string,
                               devlen);
			        }
		        }
            }
        }

		return;
}



