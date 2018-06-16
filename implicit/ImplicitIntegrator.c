/**************************************************************************
 *
 *      Module:      ImplicitIntegrator.c
 *      Description: Used with the SubcycleIntegrator function.
 *                   Implements a numerical integrator using
 *                   the implicitly-solved trapezoid integration method.
 *					 This couples the trapezoid method with the Newton-Raphson
 *					 method for nonlinear system solution. The difference
 *					 between this function and ImplicitIntegrator is that it 
 *					 only time integrates the set of forces defined by reqType. 
 *
 ***************************************************************************/
#ifdef _IMPLICIT

#include "Home.h"
#include "sys/stat.h"
#include "sys/types.h"
#include "Implicit.h"
#include "cs.h"

#ifdef _LIS
#include "lis.h"
#endif

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

#define DIR_TIMESTEP_ERROR "timestep_error"
#define DIR_NODAL_TIMESTEP "timestep"

void   SortNodes(Home_t *home, real8 *maxsep);

/*------------------------------------------------------------------------
 *
 *      Function:    UnpackQMatrices
 *      Description: This function simply unpacks the Q matrix for a node
 *      and calculates its inverse (transpose), and returns the two as
 *      the matrices Q and Qinv.
 *
 *-----------------------------------------------------------------------*/

void UnpackQMatrices(Node_t *node, real8 Q[3][3], real8 Qinv[3][3])
{
		Q[0][0] = node->Q11;
		Q[0][1] = node->Q12;
		Q[0][2] = node->Q13;
		Q[1][0] = node->Q21;
		Q[1][1] = node->Q22;
		Q[1][2] = node->Q23;
		Q[2][0] = node->Q31;
		Q[2][1] = node->Q32;
		Q[2][2] = node->Q33;

		Matrix33Transpose(Q, Qinv);
}

/*------------------------------------------------------------------------
 *
 *      Function:    ConstrainJacobian
 *      Description: Function for constraining partial Jacobians by zeroing
 *		out illegal elements and rotating them into the necessary coordinate
 *		systems.
 *
 *
 *-----------------------------------------------------------------------*/
void ConstrainJacobian(real8 Qseg[3][3], real8 Qinvseg[3][3], real8 Q1[3][3],
						real8 Qinv1[3][3], real8 Q2[3][3], real8 Qinv2[3][3],
						real8 J0[3][3], int ndof1, int ndof2, int locs[4],
						real8 J01a[3][3], real8 J01b[3][3], real8 J02a[3][3],
						real8 J02b[3][3])
{
		real8	tmpMat[3][3], tmpMat2[3][3], tmpMat3[3][3];
		int		i, j;

		//1-1 interaction (top left)
		if (locs[0]) {
			Matrix33Mult33(Qinv1, Qseg, tmpMat);
			Matrix33Mult33(Qinvseg, Q1, tmpMat2);
			Matrix33Mult33(tmpMat, J0, tmpMat3);
			Matrix33Mult33(tmpMat3, tmpMat2, J01a);
			if (ndof1==1) {
								J01a[0][1] = 0; J01a[0][2] = 0;
				J01a[1][0] = 0; J01a[1][1] = 0; J01a[1][2] = 0;
				J01a[2][0] = 0; J01a[2][1] = 0; J01a[2][2] = 0;
			} else if (ndof1==2) {
				//Do nothing
				/*
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						J01a[i][j] = J0[i][j];
					}
				} */
			} else {
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						J01a[i][j] = 0;
					}
				}
			}
		} else {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					J01a[i][j] = 0;
				}
			}
		}
							
		//1-2 interaction (top right)
		if (locs[1]) {
			Matrix33Mult33(Qinv1, Qseg, tmpMat);
			Matrix33Mult33(Qinvseg, Q2, tmpMat2);
			Matrix33Mult33(tmpMat, J0, tmpMat3);
			Matrix33Mult33(tmpMat3, tmpMat2, J01b);
			if (ndof1==1 && ndof2==1) {
								J01b[0][1] = 0; J01b[0][2] = 0;
				J01b[1][0] = 0; J01b[1][1] = 0; J01b[1][2] = 0;
				J01b[2][0] = 0; J01b[2][1] = 0; J01b[2][2] = 0;
			} else if (ndof1==1 && ndof2==2) {
				//Matrix33Mult33(Qinv1, Qseg, tmpMat);
				//Matrix33Mult33(tmpMat, J0, J01b);

				J01b[1][0] = 0; J01b[1][1] = 0;	J01b[1][2] = 0;
				J01b[2][0] = 0; J01b[2][1] = 0; J01b[2][2] = 0;
			} else if (ndof1==2 && ndof2==1) {
				//Matrix33Mult33(Qinvseg, Q2, tmpMat);
				//Matrix33Mult33(J0, tmpMat, J01b);
								J01b[0][1] = 0; J01b[0][2] = 0;
								J01b[1][1] = 0; J01b[1][2] = 0;
								J01b[2][1] = 0; J01b[2][2] = 0;
			} else if (ndof1==2 && ndof2==2) {
				//Do nothing
				/*
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						J01b[i][j] = J0[i][j];
					}
				}*/
			} else {
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						J01b[i][j] = 0;
					}
				}
			}
		} else {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					J01b[i][j] = 0;
				}
			}
		}
					
		//2-2 interaction (bottom right)
		if (locs[2]) {
			Matrix33Mult33(Qinv2, Qseg, tmpMat);
			Matrix33Mult33(Qinvseg, Q2, tmpMat2);
			Matrix33Mult33(tmpMat, J0, tmpMat3);
			Matrix33Mult33(tmpMat3, tmpMat2, J02b);
			if (ndof2==1) {
								J02b[0][1] = 0; J02b[0][2] = 0;
				J02b[1][0] = 0; J02b[1][1] = 0; J02b[1][2] = 0;
				J02b[2][0] = 0; J02b[2][1] = 0; J02b[2][2] = 0;
			} else if (ndof2==2) {
				//Do nothing
				/*
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						J02b[i][j] = J0[i][j];
					}
				} */
			} else {
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						J02b[i][j] = 0;
					}
				}
			}
		} else {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					J02b[i][j] = 0;
				}
			}
		}			

		//2-1 interaction (bottom left)
		if (locs[3]) {
			Matrix33Mult33(Qinv2, Qseg, tmpMat);
			Matrix33Mult33(Qinvseg, Q1, tmpMat2);
			Matrix33Mult33(tmpMat, J0, tmpMat3);
			Matrix33Mult33(tmpMat3, tmpMat2, J02a);
			if (ndof1==1 && ndof2==1) {
								J02a[0][1] = 0; J02a[0][2] = 0;
				J02a[1][0] = 0; J02a[1][1] = 0; J02a[1][2] = 0;
				J02a[2][0] = 0; J02a[2][1] = 0; J02a[2][2] = 0;
			} else if (ndof1==2 && ndof2==1) {
				//Matrix33Mult33(Qinv2, Qseg, tmpMat);
				//Matrix33Mult33(tmpMat, J0, J02a);

				J02a[1][0] = 0; J02a[1][1] = 0;	J02a[1][2] = 0;
				J02a[2][0] = 0; J02a[2][1] = 0; J02a[2][2] = 0;
			} else if (ndof1==1 && ndof2==2) {
				//Matrix33Mult33(Qinvseg, Q1, tmpMat);
				//Matrix33Mult33(J0, tmpMat, J02a);
								J02a[0][1] = 0; J02a[0][2] = 0;
								J02a[1][1] = 0; J02a[1][2] = 0;
								J02a[2][1] = 0; J02a[2][2] = 0;
			} else if (ndof1==2 && ndof2==2) {
				//Do nothing
				/*
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						J02a[i][j] = J0[i][j];
					}
				}*/
			} else {
				for (i=0; i<3; i++) {
					for (j=0; j<3; j++) {
						J02a[i][j] = 0;
					}
				}
			}
		} else {
			for (i=0; i<3; i++) {
				for (j=0; j<3; j++) {
					J02a[i][j] = 0;
				}
			}
		}

		return;
}	

/*------------------------------------------------------------------------
 *
 *      Function:    ConstructEnhanceList
 *      Description: Construct an enhancement list for adding elastic
 *      interaction terms to the Jacobian. This list has the following
 *      structure:
 *      -rows numbers correspond to node indices
 *      -first column is number of segments to add enhancement for each node
 *      -columns in groups of three correspond to node that forms segment
 *       with main node and the nodes of the segment within the enhancement
 *       radius of it
 *
 *      Ex: If the node with index 10 has is on a segment with node 62
 *      within renh of the segment with nodes 14 and 81, that entry
 *      would look like:
 *      row 11 (c index 10)---> ... 62 14 81 ...
 *
 *		We only store the information for a pair of segments in the row
 *		of the lowest node with J update flagged, so there is no
 *		redundancy.
 *
 *      We use the same "cell2" structure that is used for handling
 *      collisions to enable order N scaling of distance calculations.
 *      This assumes that the collision search radius with always be
 *      larger than the enhancement radius.
 *
 *      Much of this code is copied from PredictiveCollision().
 *
 *
 *-----------------------------------------------------------------------*/

void ConstructEnhanceList(Home_t *home, Node_t ***enhancelist, int *esize)//enhancelist[][4])//[][MAX_ENHANCE_NUM])
{
		int     eloc, i, j, k, q, arm12, arm21, arm34, arm43;
		int     splitSeg1, splitSeg2;
		int     cell2Index, nbrCell2Index, nextIndex;
		int     cell2X, cell2Y, cell2Z, cx, cy, cz;
		int		N, ind, arm, tag1, tag2, tag3, tag4, tagtmp;
		int		flag2, flag3, flag4;
		int		foundmatch, same13, same14, same23, same24;
		real8	renh;
		real8	tmp, dist2, dist2_1, dist2_2, ddist2dt, L1, L2;
		real8	x1, x2, x3, x4;
		real8	y1, y2, y3, y4;
		real8	z1, z2, z3, z4;
		real8	vec1[3], vec2[3];
		real8	p1[3], p2[3], p3[3], p4[3];

		//int		enhancelist[1][1];

		Node_t  *node1, *node2, *node3, *node4, *tmpNbr, *nodetmp;
		Param_t *param;

		param = home->param;

		N = home->newNodeKeyPtr;

		//NEED TO MAKE THIS AN INPUT CONTROL PARAMETER
		renh = param->renh;
		//real8 	renh = 1;
		/////////

/*
 * 		If the enhancement radius is zero or elastic interactions are
 * 		not enabled, we are done.
 */
		if (renh==0 || !param->elasticinteraction) {
			return;
		}

/*
 * 		Need to sort nodes again because nodes that have undergone collisions
 * 		and topological changes will not have cell2 membership.
 */
		SortNodes(home,&renh);

/*
 *      Start looping through native nodes looking for segments to collide...
 */
        for (tag1 = 0; tag1 < home->newNodeKeyPtr; tag1++) {

            if ((node1 = home->nodeKeys[tag1]) == (Node_t *)NULL) continue;

 /*
  * 		If the main node isn't flagged for update, skip it.
  */
			if ((node1->flags & UPDATE_NODE_J) == 0) {
				continue;
			}

/*
 * 			Add the self terms for every arm of node1
 * 			if node1 is the lowest indexed node on the arm with J
 * 			update flagged. Note we don't need
 * 			to test for proximity here - all self terms are included
 * 			for any nonzero enhancement radius.
 */
            for (arm12 = 0; arm12 < node1->numNbrs; arm12++) {

                node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);

                if (node2 == (Node_t *)NULL) continue;

                flag2 = (node2->flags & UPDATE_NODE_J) != 0;

				if ((OrderNodes(node1, node2) >= 0) || !flag2) {
					continue;
				}
				enhancelist[*esize][0] = node1;
				enhancelist[*esize][1] = node2;
				enhancelist[*esize][2] = node1;
				enhancelist[*esize][3] = node2;
				(*esize)++;
            }

/*
 *          Loop through all cell2s neighboring the node.  Only
 *          nodes in these neighboring cell2s are candidates for
 *          elastic enhancement.
 */
            cell2Index = node1->cell2Idx;
            if (cell2Index < 0) {
                continue;
            }

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

                    flag3 = (node3->flags & UPDATE_NODE_J) != 0;

/*
 *                  Loop over all arms of node1.
 */
                    for (arm12 = 0; arm12 < node1->numNbrs; arm12++) {

                        node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);

                        if (node2 == (Node_t *)NULL) continue;

                        flag2 = (node2->flags & UPDATE_NODE_J) != 0;

                        if ((OrderNodes(node1, node2) >= 0) || !flag2) {
                            continue;
                        }

/*
 *                      Loop over all arms of node3.
 */
                        for (arm34 = 0; arm34 < node3->numNbrs; arm34++) {

                            node4 = GetNodeFromTag(home, node3->nbrTag[arm34]);

                            if (node4 == (Node_t *)NULL) continue;

/*
 *                      	Ensures the node with the lower tag is the node3
 */
							if (OrderNodes(node3, node4) >= 0) {
								continue;
							}

/*
 *                      	Make sure this isn't a self term (handled above),
 *                      	and that we only do this segment pair once.
 */
							if ((node1 == node3) && (node2 == node4)) {
								continue;
							}

							if ((OrderNodes(node3, node1) < 0) || !flag3) {
								continue;
							}

							if ((OrderNodes(node4, node2) < 0) &&
								(node1 == node3)) {
								continue;
							}

                            x1 = node1->x; y1 = node1->y; z1 = node1->z;
                            x2 = node2->x; y2 = node2->y; z2 = node2->z;
                            x3 = node3->x; y3 = node3->y; z3 = node3->z;
                            x4 = node4->x; y4 = node4->y; z4 = node4->z;

                            PBCPOSITION(param, x1, y1, z1, &x2, &y2, &z2);
                            PBCPOSITION(param, x1, y1, z1, &x3, &y3, &z3);
                            PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

                            p1[0] = x1;  p1[1] = y1;  p1[2] = z1;
                            p2[0] = x2;  p2[1] = y2;  p2[2] = z2;
                            p3[0] = x3;  p3[1] = y3;  p3[2] = z3;
                            p4[0] = x4;  p4[1] = y4;  p4[2] = z4;

/*
 *                          It is possible to have a zero-length segment.
 *                          If encountered, avoid it.
 */
                            vec1[0] = x2 - x1;
                            vec1[1] = y2 - y1;
                            vec1[2] = z2 - z1;

                            if (DotProduct(vec1, vec1) < 1.0e-20) {
                                continue;
                            }

                            vec2[0] = x4 - x3;
                            vec2[1] = y4 - y3;
                            vec2[2] = z4 - z3;

                            if (DotProduct(vec2, vec2) < 1.0e-20) {
                                continue;
                            }
/*
 *                          Find the minimum distance between the two segments
 *                          and determine if their enhancements are necessary.
 */
							GetMinDist2(x1, y1, z1, 0, 0, 0,
										x2, y2, z2, 0, 0, 0,
										x3, y3, z3, 0, 0, 0,
										x4, y4, z4, 0, 0, 0,
										&dist2, &ddist2dt, &L1, &L2);


/*
 * 							Add information for all nodes - we already checked
 * 							above if this pair of segments has been covered.
 */
							if (sqrt(dist2)<renh) {
								enhancelist[*esize][0] = node1;
								enhancelist[*esize][1] = node2;
								enhancelist[*esize][2] = node3;
								enhancelist[*esize][3] = node4;
								(*esize)++;
							}

                        } //for (arm34)
                    } //for (arm12)
                } //while (nextnode)
              } //for (cz)
             } //for (cy)
            } //for (cx)
        } //for (tag1=main node)
        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    CalcPerturbedForces
 *      Description: This function is used for calculating Jacobian
 *      elastic interaction terms using finite differences. Given 4 nodes,
 *      it calculates the force change due to perturbing one of them
 *      (nodep) in its local x and y directions. The difference between
 *      these changes is then passed out, to be used by CalcJElTerms
 *      function. If two of the nodes are the same (node1=node3 and
 *      node2=node4), this is interpretted as only using the elastic self
 *      force change.
 *
 *-----------------------------------------------------------------------*/

void CalcPerturbedForces(Home_t *home, Node_t *node1, Node_t *node2, Node_t *node3,
						   Node_t *node4, Node_t *nodep, real8 delta,
						   real8 df1x[3], real8 df2x[3], real8 df3x[3], real8 df4x[3],
						   real8 df1y[3], real8 df2y[3], real8 df3y[3], real8 df4y[3])
{

		int 	selfterm, arm;
		real8	MU, NU, a, bx, by, bz;
		real8	r[3], r_loc[3], r_tmp[3];
		real8	f1p[3], f2p[3], f3p[3], f4p[3];
		real8	f1m[3], f2m[3], f3m[3], f4m[3];
		real8	Qp[3][3], Qinvp[3][3];
		Param_t *param;

		param = home->param;

/*
 * 		If the perturbed node has no dofs, return zeros.
 */
		if (nodep->ndof==0) {
			df1x[0] = 0; df1x[1] = 0; df1x[2] = 0;
			df2x[0] = 0; df2x[1] = 0; df2x[2] = 0;
			df3x[0] = 0; df3x[1] = 0; df3x[2] = 0;
			df4x[0] = 0; df4x[1] = 0; df4x[2] = 0;

			df1y[0] = 0; df1y[1] = 0; df1y[2] = 0;
			df2y[0] = 0; df2y[1] = 0; df2y[2] = 0;
			df3y[0] = 0; df3y[1] = 0; df3y[2] = 0;
			df4y[0] = 0; df4y[1] = 0; df4y[2] = 0;
			return;
		}

/*
 * 		Unpack the necessary Q matrices.
 */
		UnpackQMatrices(nodep, Qp, Qinvp);

/*
 * 		Check if this is a self force term. If it is, unpack a few
 * 		things.
 */
		if ((node1->myTag.domainID == node3->myTag.domainID)&&
				(node1->myTag.index == node3->myTag.index) &&
				(node2->myTag.domainID == node4->myTag.domainID)&&
				(node2->myTag.index == node4->myTag.index)) {
			selfterm = 1;
			MU = param->shearModulus;
			NU = param->pois;
			a = param->rc;
			arm = GetArmID (home, node1, node2);
			bx = node1->burgX[arm];
			by = node1->burgY[arm];
			bz = node1->burgZ[arm];
		} else {
			selfterm = 0;
		}

/*
 * 		Calculate the local coordinates of the node to be perturbed.
 */
		r[0] = nodep->x;
		r[1] = nodep->y;
		r[2] = nodep->z;

		Matrix33Vector3Multiply(Qinvp, r, r_loc);


/*
 * 		Perturb the node in its local x-direction, rotate into global
 * 		coordinates, and set as new nodal position. Also, reset r_loc
 * 		vector.
 */

		r_loc[0] = r_loc[0] + delta;
		Matrix33Vector3Multiply(Qp, r_loc, r_tmp);
		r_loc[0] = r_loc[0] - delta;

		nodep->x = r_tmp[0];
		nodep->y = r_tmp[1];
		nodep->z = r_tmp[2];

/*
 * 		Calculate nodal forces in the perturbed configuration.
 */
		if (selfterm) {
			SelfForce(0, MU, NU, bx, by, bz, node1->x, node1->y, node1->z,
							node2->x, node2->y, node2->z, a, 0, f1p, f2p);
			f3p[0] = 0; f3p[1] = 0; f3p[2] = 0;
			f4p[0] = 0; f4p[1] = 0; f4p[2] = 0;
		} else {
			ComputeForces(home, node1, node2, node3, node4, f1p, f2p, f3p, f4p);
		}

/*
 * 		If the perturbed node has no y dof, return zeros for dfy.
 */
		if (nodep->dofy==0) {
			df1y[0] = 0; df1y[1] = 0; df1y[2] = 0;
			df2y[0] = 0; df2y[1] = 0; df2y[2] = 0;
			df3y[0] = 0; df3y[1] = 0; df3y[2] = 0;
			df4y[0] = 0; df4y[1] = 0; df4y[2] = 0;
			
			nodep->x = r[0];
			nodep->y = r[1];
			nodep->z = r[2];
			return;
		}

/*
 * 		Perturb in the other direction and repeat.
 */
		r_loc[0] = r_loc[0] - delta;
		Matrix33Vector3Multiply(Qp, r_loc, r_tmp);
		r_loc[0] = r_loc[0] + delta;

		nodep->x = r_tmp[0];
		nodep->y = r_tmp[1];
		nodep->z = r_tmp[2];

		if (selfterm) {
			SelfForce(0, MU, NU, bx, by, bz, node1->x, node1->y, node1->z,
							node2->x, node2->y, node2->z, a, 0, f1m, f2m);
			f3m[0] = 0; f3m[1] = 0; f3m[2] = 0;
			f4m[0] = 0; f4m[1] = 0; f4m[2] = 0;
		} else {
			ComputeForces(home, node1, node2, node3, node4, f1m, f2m, f3m, f4m);
		}

/*
 * 		Calculate delta-f for the finite difference calculation.
 */
		df1x[0] = f1p[0]-f1m[0]; df1x[1] = f1p[1]-f1m[1]; df1x[2] = f1p[2]-f1m[2];
		df2x[0] = f2p[0]-f2m[0]; df2x[1] = f2p[1]-f2m[1]; df2x[2] = f2p[2]-f2m[2];
		df3x[0] = f3p[0]-f3m[0]; df3x[1] = f3p[1]-f3m[1]; df3x[2] = f3p[2]-f3m[2];
		df4x[0] = f4p[0]-f4m[0]; df4x[1] = f4p[1]-f4m[1]; df4x[2] = f4p[2]-f4m[2];

/*
 * 		Repeat everything for perturbations in the local y-direction.
 */
		r_loc[1] = r_loc[1] + delta;
		Matrix33Vector3Multiply(Qp, r_loc, r_tmp);
		r_loc[1] = r_loc[1] - delta;

		nodep->x = r_tmp[0];
		nodep->y = r_tmp[1];
		nodep->z = r_tmp[2];

		if (selfterm) {
			SelfForce(0, MU, NU, bx, by, bz, node1->x, node1->y, node1->z,
							node2->x, node2->y, node2->z, a, 0, f1p, f2p);
			f3p[0] = 0; f3p[1] = 0; f3p[2] = 0;
			f4p[0] = 0; f4p[1] = 0; f4p[2] = 0;
		} else {
			ComputeForces(home, node1, node2, node3, node4, f1p, f2p, f3p, f4p);
		}

		r_loc[1] = r_loc[1] - delta;
		Matrix33Vector3Multiply(Qp, r_loc, r_tmp);
		r_loc[1] = r_loc[1] + delta;

		nodep->x = r_tmp[0];
		nodep->y = r_tmp[1];
		nodep->z = r_tmp[2];

		if (selfterm) {
			SelfForce(0, MU, NU, bx, by, bz, node1->x, node1->y, node1->z,
							node2->x, node2->y, node2->z, a, 0, f1m, f2m);
			f3m[0] = 0; f3m[1] = 0; f3m[2] = 0;
			f4m[0] = 0; f4m[1] = 0; f4m[2] = 0;
		} else {
			ComputeForces(home, node1, node2, node3, node4, f1m, f2m, f3m, f4m);
		}

		df1y[0] = f1p[0]-f1m[0]; df1y[1] = f1p[1]-f1m[1]; df1y[2] = f1p[2]-f1m[2];
		df2y[0] = f2p[0]-f2m[0]; df2y[1] = f2p[1]-f2m[1]; df2y[2] = f2p[2]-f2m[2];
		df3y[0] = f3p[0]-f3m[0]; df3y[1] = f3p[1]-f3m[1]; df3y[2] = f3p[2]-f3m[2];
		df4y[0] = f4p[0]-f4m[0]; df4y[1] = f4p[1]-f4m[1]; df4y[2] = f4p[2]-f4m[2];

/*
 * 		Restore the perturbed node to its original position.
 */
		nodep->x = r[0];
		nodep->y = r[1];
		nodep->z = r[2];

}

/*------------------------------------------------------------------------
 *
 *      Function:    CalcJElterms
 *      Description: Calculate a sub-Jacobian of elastic interaction terms.
 *      We assume a linear, glide-only mobility law is being used.
 *      node = node whose velocity is changing (row)
 *      dfx and dfy = vectors of force differences due to perturbing a node.
 *
 *-----------------------------------------------------------------------*/

void CalcJElterms(Node_t *node, real8 dfx[3], real8 dfy[3], real8 delta,
					real8 B, real8 Lsum, real8 subJ[3][3])
{
		real8	Q[3][3], Qinv[3][3];
		real8	dfx_loc[3], dfy_loc[3];
		real8 	factor;
/*
 * 		Unpack the necessary Q matrices.
 */
		UnpackQMatrices(node, Q, Qinv);

/*
 * 		Transform force difference into local coordinate system of
 * 		response node.
 */
		Matrix33Vector3Multiply(Qinv, dfx, dfx_loc);
		Matrix33Vector3Multiply(Qinv, dfy, dfy_loc);

/*
 * 		Calculate the sub-Jacobian.
 */
		factor = 1/(2*delta*B*(Lsum/2));
		subJ[0][0] = factor*dfx_loc[0];
		subJ[0][1] = factor*dfx_loc[1];
		subJ[0][2] = factor*dfx_loc[2];
		subJ[1][0] = factor*dfy_loc[0];
		subJ[1][1] = factor*dfy_loc[1];
		subJ[1][2] = factor*dfy_loc[2];
		subJ[2][0] = 0;
		subJ[2][1] = 0;
		subJ[2][2] = 0;
}


/*------------------------------------------------------------------------
 *
 *      Function:    AddJterms
 *      Description: Add new Jacobian terms to the Jacobian arrays.
 *      We first check to see if a J element exists, and if it does we
 *      add to it, otherwise we create a new element.
 *      node1 = node whose velocity is changing
 *      node2 = node being perturbed
 *
 *-----------------------------------------------------------------------*/

void AddJterms(Param_t *param, Node_t *node1, Node_t *node2, real8 subJ[3][3], cs *Jmob)
{
		int		indrow[3], indcol[3];
		int		i, j, k, foundit, flag1, flag2;
		real8	tol = 1e-8;
		int 	err;

		indrow[0] = node1->dofx;
		indrow[1] = node1->dofy;
		indrow[2] = node1->dofz;

		indcol[0] = node2->dofx;
		indcol[1] = node2->dofy;
		indcol[2] = node2->dofz;

////////////////////////////////////////
/*
 *		If we are using an FCC-type mobility law, any node on a non-(111)
 *		glide plane can only move along its line - this corresponds to a zero-
 *		energy mode, and hence the associated Jacobian terms are zero. So, we
 *		need to check for this and prevent adding Jterms for any such node.
 */
/*
		if (param->materialType==MAT_TYPE_FCC) {
			for (i=0; i<node1->numNbrs; i++) {
				if (fabs(fabs(node1->nx[i])-fabs(node1->ny[i]))>tol ||
						fabs(fabs(node1->ny[i])-fabs(node1->nz[i]))>tol) {
					indrow[0] = 0;
					indrow[1] = 0;
					indrow[2] = 0;
					break;
				}
			}

			for (i=0; i<node2->numNbrs; i++) {
				if (fabs(fabs(node2->nx[i])-fabs(node2->ny[i]))>tol ||
						fabs(fabs(node2->ny[i])-fabs(node2->nz[i]))>tol) {
					indcol[0] = 0;
					indcol[1] = 0;
					indcol[2] = 0;
					break;
				}
			}
		}
*/
////////////////////////////////////////

		flag1 = ((node1->flags & UPDATE_NODE_J)!=0);
		flag2 = ((node2->flags & UPDATE_NODE_J)!=0);

		for (i=0; i<3; i++) {
			for (j=0; j<3; j++) {
				//If the element is close to zero, omit it.
				if (fabs(subJ[i][j])<tol) {
					continue;
				}
				foundit = 0;
				//If one of the indices is zero, does not correspond to a
				//dof - don't add these contributions. Also, at least one node
				//must be flagged.
				if ((indrow[i]!=0 && indcol[j]!=0)
						&& (flag1 || flag2)){
					cs_entry(Jmob, (csi)indrow[i]-1, (csi)indcol[j]-1,(double)subJ[i][j]);
//printf("this element is %i %i %e\n",Jmob->p[Jmob->nz-1],Jmob->i[Jmob->nz-1],Jmob->x[Jmob->nz-1]);
/*
					for (k=0; k<*Jsize; k++) {
						if ((Jrow[k]==indrow[i]) &&
								(Jcol[k]==indcol[j])) {
							Jval[k] = Jval[k] + subJ[i][j];
							foundit = 1;
							break;
						}
					}
					if (!foundit) {
						Jrow[*Jsize] = indrow[i];
						Jcol[*Jsize] = indcol[j];
						Jval[*Jsize] = subJ[i][j];
						(*Jsize)++;
					}
*/
				}
			}
		}
}


/*------------------------------------------------------------------------
 *
 *      Function:    EnhanceJacobian
 *      Description: Given four nodes, calculate and add all elastic
 *      interaction terms associated with them to the Jacobian arrays.
 *
 *-----------------------------------------------------------------------*/

void EnhanceJacobian(Home_t *home, Node_t *node1, Node_t *node2, Node_t *node3,
						Node_t *node4, real8 B, real8 Lsum1, real8 Lsum2,
						real8 Lsum3, real8 Lsum4, cs *Jmob)
{
		int		index1, index2, index3, index4, selfterm;

		real8	delta;
		real8	df1x[3], df2x[3], df3x[3], df4x[3];
		real8	df1y[3], df2y[3], df3y[3], df4y[3];
		real8	subJ[3][3];
		Param_t *param;

		param = home->param;

/*
 * 		Set the perturbation distance for elastic interaction calculation.
 */
		delta = 1;

//		index1 = node1->myTag.index;
//		index2 = node2->myTag.index;
//		index3 = node3->myTag.index;
//		index4 = node4->myTag.index;

/*
 * 		Determine if this is a self interaction.
 */
		if (node1==node3 && node2==node4) {
			selfterm = 1;
		} else {
			selfterm = 0;
		}

/*
 * 		Calculate the terms for each nodal perturbation, and add them to
 * 		the Jacobian arrays.
 */
		//Perturb node 1
		CalcPerturbedForces(home, node1, node2, node3, node4, node1,
					delta, df1x, df2x, df3x, df4x, df1y, df2y, df3y, df4y);

		CalcJElterms(node1, df1x, df1y, delta, B, Lsum1, subJ);
		AddJterms(param, node1, node1, subJ, Jmob);

		CalcJElterms(node2, df2x, df2y, delta, B, Lsum2, subJ);
		AddJterms(param, node2, node1, subJ, Jmob);

		if (!selfterm) {
			CalcJElterms(node3, df3x, df3y, delta, B, Lsum3, subJ);
			AddJterms(param, node3, node1, subJ, Jmob);

			CalcJElterms(node4, df4x, df4y, delta, B, Lsum4, subJ);
			AddJterms(param, node4, node1, subJ, Jmob);
		}

		//Perturb node 2
		CalcPerturbedForces(home, node1, node2, node3, node4, node2,
					delta, df1x, df2x, df3x, df4x, df1y, df2y, df3y, df4y);
		CalcJElterms(node1, df1x, df1y, delta, B, Lsum1, subJ);
		AddJterms(param, node1, node2, subJ, Jmob);

		CalcJElterms(node2, df2x, df2y, delta, B, Lsum2, subJ);
		AddJterms(param, node2, node2, subJ, Jmob);

		if (!selfterm) {
			CalcJElterms(node3, df3x, df3y, delta, B, Lsum3, subJ);
			AddJterms(param, node3, node2, subJ, Jmob);

			CalcJElterms(node4, df4x, df4y, delta, B, Lsum4, subJ);
			AddJterms(param, node4, node2, subJ, Jmob);
		}


/*
 * 		If this is a self interaction, we don't need to perturb the same
 * 		nodes again. Also, if these two segments share a node, we don't
 * 		need to perturb the common node again.
 */
		if (!selfterm) {
			if (node1!=node3 && node2!=node3) {
				//Perturb node 3
				CalcPerturbedForces(home, node1, node2, node3, node4, node3,
							delta, df1x, df2x, df3x, df4x, df1y, df2y, df3y, df4y);

				CalcJElterms(node1, df1x, df1y, delta, B, Lsum1, subJ);
				AddJterms(param, node1, node3, subJ, Jmob);

				CalcJElterms(node2, df2x, df2y, delta, B, Lsum2, subJ);
				AddJterms(param, node2, node3, subJ, Jmob);

				CalcJElterms(node3, df3x, df3y, delta, B, Lsum3, subJ);
				AddJterms(param, node3, node3, subJ, Jmob);

				//Check that we haven't already covered this case above.
				if (node1!=node4 && node2!=node4) {
					CalcJElterms(node4, df4x, df4y, delta, B, Lsum4, subJ);
					AddJterms(param, node4, node3, subJ, Jmob);
				}
			}

			if (node1!=node4 && node2!=node4) {
				//Perturb node 4
				CalcPerturbedForces(home, node1, node2, node3, node4, node4,
							delta, df1x, df2x, df3x, df4x, df1y, df2y, df3y, df4y);

				CalcJElterms(node1, df1x, df1y, delta, B, Lsum1, subJ);
				AddJterms(param, node1, node4, subJ, Jmob);

				CalcJElterms(node2, df2x, df2y, delta, B, Lsum2, subJ);
				AddJterms(param, node2, node4, subJ, Jmob);

				//Check that we haven't already covered this case above.
				if (node1!=node3 && node2!=node3) {
					CalcJElterms(node3, df3x, df3y, delta, B, Lsum3, subJ);
					AddJterms(param, node3, node4, subJ, Jmob);
				}

				CalcJElterms(node4, df4x, df4y, delta, B, Lsum4, subJ);
				AddJterms(param, node4, node4, subJ, Jmob);
			}
		}
}


/*------------------------------------------------------------------------
 *
 *      Function:    AssignCoordSys
 *      Description: Assigns a coodinate system to each node
 *      for the implicit solver. Node coordinate systems are assigned as follows:
 *      	1) If node arms only have one glide plane, node is assignment
 *      	coordinate system of the arm it shares with neighbor of lowest
 *      	tag number.
 *      	2) If node arms have more than one glide plane, coordinate system
 *      	is set by glide constraints.
 *
 *      The dof_list is also calculated in this function. It stores degree
 *      of freedom number assignments for each nodal degree of freedom.
 *
 *-----------------------------------------------------------------------*/

void AssignCoordSys(Home_t *home, int *Ndof, int group)
{

	real8	 nxtmp, nytmp, nztmp, absdot, nxdiff, nydiff, nzdiff, tmp;
	real8	 dx, dy, dz, linex, liney, linez;
	int		 N, i, j, k, l, numplanes, ndofi, nbrindex, low_nbr, low_index;
	int		 dofnum, ind, cycle, on_non111, inSegGroup;
	real8	 planes[100]; //arbitrarily large
	int		 thisdof, olddofs[3], Jrowold[10000], Jcolold[10000];
	real8	 xdir[3], ydir[3], zdir[3], tmpvec1[3], tmpvec2[3];
	real8 	 tol = 1e-5;
	Node_t	 *node1 , *node2 , *node3 , *node4 , *node, *nbr;
	Param_t	 *param;
	Subcyc_t *subcyc;

	Implicit_t *implicit;
	implicit = home->implicit;
	param    = home->param;
	subcyc   = home->subcyc;

	N = home->newNodeKeyPtr;

/*
 *	Loop through all nodes. If a node is on the Jacobian update list,
 *	calculate its transformation matrix.
 */
	for (i=0; i < N; i++) {
		if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
		
		if (subcyc->SegListG0_cnt + subcyc->SegSegListG0_cnt == 0) {
			inSegGroup = 1;
		} else {
			inSegGroup = 0;
			
			for (j = 0; j < subcyc->SegListG0_cnt + subcyc->SegSegListG0_cnt ; j++) {
				if (j < subcyc->SegListG0_cnt) {
					if (subcyc->SegListG0[j].flag == 0) continue;
					
					node1 = subcyc->SegListG0[j].seg->node1;
					node2 = subcyc->SegListG0[j].seg->node2;
					
					if (node==node1  || node==node2) {
						inSegGroup = 1;
						break;
					}
				} else {
					k = j - subcyc->SegListG0_cnt;
					if (subcyc->SegSegListG0[k].flag == 0) continue;
				
					node1 = subcyc->SegSegListG0[k].seg1->node1;
					node2 = subcyc->SegSegListG0[k].seg1->node2;
					node3 = subcyc->SegSegListG0[k].seg2->node1;
					node4 = subcyc->SegSegListG0[k].seg2->node2;
					
					if (node==node1  || node==node2 || node==node3  || node==node4) {
						inSegGroup = 1;
						break;
					}
				}
			}
		}
		
/*
 *		If node is not in the group, skip it.
 */
		if (inSegGroup==0 || (group >= 0 && node->subgroup != group)) {
			node->Q11 = 1;   node->Q21 = 0;   node->Q31 = 0;
			node->Q12 = 0;   node->Q22 = 1;   node->Q32 = 0;
			node->Q13 = 0;   node->Q23 = 0;   node->Q33 = 1;
			node->ndof = 0;
			continue;
		}

		if ((node->flags & UPDATE_NODE_J) == 0) continue;
/*
 * 		Determine how many unique glide planes the node is on.
 */
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
 *		If we are using an FCC-type mobility law, any node on a non-(111)
 *		glide plane can only move along its line. First we check to see if the
 *		node is on any non-(111) planes. If it is, we count the number of 
 *		unique line directions its arms have. If more than one, the node cannot
 *		move. If it only has one line direction, we finally test to make sure
 *		that direction lies in all glide planes.  
 */
		on_non111 = 0;		
		if (param->materialType==MAT_TYPE_FCC) {
			for (k=0; k<node->numNbrs; k++) {
				if (fabs(fabs(node->nx[k])-fabs(node->ny[k]))>tol ||
						fabs(fabs(node->ny[k])-fabs(node->nz[k]))>tol) {
					on_non111 = 1;
					break;
				}
			}
			if (on_non111) {
				//Test for more than one unique line direction.
				for (k=0; k<node->numNbrs; k++) {
					nbr = GetNeighborNode(home,node,k);					
					dx = node->x-nbr->x;
					dy = node->y-nbr->y;
					dz = node->z-nbr->z;
					ZImage(param,&dx,&dy,&dz);
					Normalize(&dx,&dy,&dz);
					if (k==0) {
						linex = dx; liney = dy; linez = dz;
					} else {
						absdot = fabs(dx*linex+dy*liney+dz*linez);
						if (fabs(1.0-absdot)>tol) {					
							numplanes = 3; //no dofs assigned below.
							break;
						}
					}
				}
				//If only one line direction, test it against the glide planes.
				if (numplanes<3) {
					for (k=0; k<numplanes; k++) {
						absdot = fabs(linex*planes[3*k]+
										liney*planes[3*k+1]+
											linez*planes[3*k+2]);
						if (absdot>tol) {					
							numplanes = 3; //no dofs assigned below
							break;
						}
					}
				}
			}
		}
						
			
/*
 * 		If the node is on more than two planes, it is immobile and needs no Q matrix.
 * 		If it is on two planes, it only has one degree of freedom and the x-axis
 * 		is its only dof. If it is on one plane, it has two dofs.
 */
		if (numplanes>2) {
			xdir[0] = 1; xdir[1] = 0; xdir[2] = 0;
			ydir[0] = 0; ydir[1] = 1; ydir[2] = 0;
			zdir[0] = 0; zdir[1] = 0; zdir[2] = 1;
			ndofi = 0;
		} else {
/*
 *			If this is an FCC simulation, and we found a non-(111) glide plane
 *			above, then the only degree of freedom is along the line direction.
 *			Otherwise, use plane information to assign a coordinate system. 
 */
			if (param->materialType==MAT_TYPE_FCC && on_non111) {
				xdir[0] = linex;
				xdir[1] = liney;
				xdir[2] = linez;
				NormalizeVec(xdir);
				zdir[0] = node->nx[0];
				zdir[1] = node->ny[0];
				zdir[2] = node->nz[0];
				NormalizeVec(zdir);
				NormalizedCrossVector(zdir, xdir, ydir);
				ndofi = 1;
			} else if (numplanes==2) {
				tmpvec1[0] = planes[0]; tmpvec1[1] = planes[1]; tmpvec1[2] = planes[2];
				tmpvec2[0] = planes[3]; tmpvec2[1] = planes[4]; tmpvec2[2] = planes[5];
				NormalizedCrossVector(tmpvec1, tmpvec2, xdir);
				ydir[0] = planes[0]; ydir[1] = planes[1]; ydir[2] = planes[2];
				NormalizedCrossVector(xdir, ydir, zdir);
				ndofi = 1;
			} else {			
/*
 * 				Find the neighbor (nbr) with the lowest index and use its
 * 				coordinate system.
 */
				for (k = 0; k < node->numNbrs; k++) {
					nbr = GetNeighborNode(home, node, k);
					nbrindex = nbr->myTag.index;
					if (k==0) {
						low_nbr = 0;
						low_index = nbrindex;
					} else {
						if (nbrindex<low_index) {
							low_nbr = k;
							low_index = nbrindex;
						}
					}
				}
				xdir[0] = node->burgX[low_nbr];
				xdir[1] = node->burgY[low_nbr];
				xdir[2] = node->burgZ[low_nbr];
				NormalizeVec(xdir);
				zdir[0] = node->nx[low_nbr];
				zdir[1] = node->ny[low_nbr];
				zdir[2] = node->nz[low_nbr];
				NormalizeVec(zdir);
				if (fabs(DotProduct(xdir,zdir))>1e-10) {
					printf("Warning: Burgers vector is not orthogonal to" 
							"the glide plane normal for node %i!\nb=%e %e %e"
							"n=%e %e %e\n",node->myTag.index,xdir[0],xdir[1],
							xdir[2],zdir[0],zdir[1],zdir[2]);
				}
				NormalizedCrossVector(zdir, xdir, ydir);
				ndofi = 2;
			}
		}

		node->Q11 = xdir[0];
		node->Q21 = xdir[1];
		node->Q31 = xdir[2];
		node->Q12 = ydir[0];
		node->Q22 = ydir[1];
		node->Q32 = ydir[2];
		node->Q13 = zdir[0];
		node->Q23 = zdir[1];
		node->Q33 = zdir[2];
		node->ndof = ndofi;


/*
 * 		If the nodes if frozen, set the number of dofs to zero.
 */
		if (node->constraint == PINNED_NODE) {
			node->ndof = 0;
		}
	}

/*
 * 	Update the dof numbers for each node. First store the old dof numbers,
 * 	update them, and then update the Jacobian arrays to reflect the dof changes.
 *
 */
	dofnum = 1;
	for (i=0; i<N; i++) {
		node = home->nodeKeys[i];
		if (node == (Node_t *)NULL) continue;
		olddofs[0] = node->dofx;
		olddofs[1] = node->dofy;
		olddofs[2] = node->dofz;
		for (j=0; j<3; j++) {
			if (j<node->ndof) {
				if (j==0){
					node->dofx = dofnum;
				} else if (j==1) {
					node->dofy = dofnum;
				} else {
					node->dofz = dofnum;
				}
				thisdof = dofnum;
				dofnum ++;
			} else {
				if (j==0){
					node->dofx = 0;
				} else if (j==1) {
					node->dofy = 0;
				} else {
					node->dofz = 0;
				}
				thisdof = 0;
			}
		}
	}
	

/*
 * 	Update the total number of dofs in the system.
 */
	*Ndof = dofnum-1;
}


/*------------------------------------------------------------------------
 *
 *      Function:    UpdateJacobian
 *      Description: Update the mobility law Jacobian, Jmob. Jmob is stored
 *		using the coordinate (COO) format whereby every element is defined by
 *		three vectors - one defining the row, Jrow, one defining the comlumn
 *		Jcol, and one defining the value, Jval. For details see the lis
 *		user's manual. 
 *
 *
 *-----------------------------------------------------------------------*/
void UpdateJacobian(Home_t *home, int newJ, int *changedJ,	cs *Jmob)
{		
		int		firststep, n, xonly1, xonly2, numplanes, add2, done1, NJ;
		int 	i, j, k, l, i1, i2, j1, rowlim, indrow, indcol, foundit, ndof1, ndof2;
		int		doneLoc, p;
		int		numdiff, numplanes1, numplanes2, locs[4], ind[6], nbrindex;
		int		maxenhnum, numer, denom, eint, esize;
		real8	Lsum, Lsum1, Lsum2, bmag2, delta;
		real8	xtmp, ytmp, ztmp;
		real8	xdiff, ydiff, zdiff;
		real8	nxtmp, nytmp, nztmp;
		real8	nxdiff, nydiff, nzdiff;
		real8	xdir1[3], ydir1[3], zdir1[3];
		real8	xdir2[3], ydir2[3], zdir2[3];
		real8	xdirseg[3], ydirseg[3], zdirseg[3];	
		real8	Qseg[3][3], Qinvseg[3][3];
		real8	Q1[3][3], Qinv1[3][3];	
		real8	Q2[3][3], Qinv2[3][3];
		real8	planes1[30], planes2[30];
		real8	vectmp[3];
		real8	stressMat[3][3], stressMatloc[3][3];
		real8	tmpMat[3][3], tmpMat2[3][3], J0[3][3];
		real8	J01a[3][3], J01b[3][3], J02a[3][3], J02b[3][3];
		real8	JL[6][6], Japp[6][6], Jcore[6][6], Jseg[6][6];
		real8	r1_loc[3], r2_loc[3], dr[3], b_loc[3], dr_norm[3]; 
		real8	f1_loc[3], f2_loc[3], fapp[3];
		real8	L, coreFactor, s, B, Ecore;

		real8	Lsum3, Lsum4, subJ[3][3];
		Node_t	*node, *nbr, *nbrtmp, *node1, *node2, *node3, *node4;
		Param_t *param;
		
        Implicit_t *implicit;
        implicit = home->implicit;

		real8 	tol = 1e-10;
		int  	N = home->newNodeKeyPtr;
		param = home->param;

		real8 	*Lsumvec = malloc((N+1) * sizeof(real8));

/*
 * 		Allocate memory for the enhancelist. The maximum number of elements
 * 		in the list is n!/(2(n-2)!) + n, where n is the number of segments.
 */
		maxenhnum = (int)((real8)N+1)*((real8)N/2+1);
		Node_t 	***enhancelist = malloc(maxenhnum * sizeof(Node_t **));
		for(i=0;i<maxenhnum;i++){
			enhancelist[i] = malloc(4 * sizeof(Node_t *));
		}


/*
 *		Check to see if this is the first time step.
 */
		firststep = 0;
 

/*
 *		If we are updating the whole Jacobian, flag every node.
 */

		if (firststep || newJ) { 
			for (i=0; i < N; i++) {
		        node = home->nodeKeys[i];
		        if (node == (Node_t *)NULL) continue;
				node->flags |= UPDATE_NODE_J;
			}
		}


/*
 * 		Construct an enhancement list for adding elastic interaction terms
 * 		to the Jacobian.
 */
		esize = 0;
		ConstructEnhanceList(home,enhancelist,&esize);

/*
 * 		Calculate the segment length sum for all nodes that will be updated.
 */
		for (i=0; i < N; i++) {
	        node = home->nodeKeys[i];
	        if (node == (Node_t *)NULL) continue;
				Lsum = 0;
				for (k = 0; k < node->numNbrs; k++) {
					nbrtmp = GetNeighborNode(home, node, k);
					xtmp = nbrtmp->x;
					ytmp = nbrtmp->y;
					ztmp = nbrtmp->z;
					PBCPOSITION(param, node->x, node->y, node->z, &xtmp, &ytmp, &ztmp);
					xdiff = node->x - xtmp;
					ydiff = node->y - ytmp;
					zdiff = node->z - ztmp;
					Lsum = Lsum + sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff);
				}
				Lsumvec[node->myTag.index] = Lsum;
		}

/*
 *		Loop through the nodes and update the Jacobian. Nodes that have
 *		undergone topological changes (remeshing or collisions/splitting)
 *		need to be updated. When changing a Jmob element, we must first check
 * 		to see if it is already defined before appending the Jmob arrays to 
 *		add its value. 
 */

/* 
 *		We will need the current stress tensor later, so unpack it here.
 */

		stressMat[0][0] = param->appliedStress[0];
		stressMat[1][1] = param->appliedStress[1];
		stressMat[2][2] = param->appliedStress[2];
		stressMat[0][1] = param->appliedStress[5];
		stressMat[0][2] = param->appliedStress[4];
		stressMat[1][0] = param->appliedStress[5];
		stressMat[1][2] = param->appliedStress[3];
		stressMat[2][0] = param->appliedStress[4];
		stressMat[2][1] = param->appliedStress[3];

/*
 *		If a node and its neighbor are flagged, we add the contributions for 
 *		that segment. If only the main node is flagged, we add only its
 *		contribution.
 */
		for (i=0; i < N; i++) {
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
			if ((node->flags & UPDATE_NODE_J) == 0) {
				continue;
			}
			for (j = 0; j < node->numNbrs; j++) {	
				nbr = GetNeighborNode(home, node, j);
/*
 *				Check to see if we have already handled this segment.
 */
				if ((OrderNodes(node, nbr) >= 0) ||
						((nbr->flags & UPDATE_NODE_J) == 0)) {
					continue;
				}

				int tag1 = node->myTag.index;
				int tag2 = nbr->myTag.index;

/*
*				Calculate the lengths of all segments attached to both
*				nodes.
*/
	    		Lsum1 = Lsumvec[node->myTag.index];
	    		Lsum2 = Lsumvec[nbr->myTag.index];

/*
 *				Calculate the local coordinate systems for both nodes
 *				and the segment.
 *				If a node only has one glide plane then the coordinate
 *				directions are as follows:
 *				x = Burgers vector direction
 *				z = glide plane normal
 *				y = cross product of z with x
 *				This is also the convention for the segment.
 *				If a node has more than one glide plane, then define
 *				the coordinate directions as follows:
 *				x = cross product of any two glide plane normals
 *				y = normal to any glide plane
 *				z = cross product of x with y
 *				Also calculate the resulting transformation matrices
 *				and store them for later use.
 */
				//Segment
				xdirseg[0] = node->burgX[j];
				xdirseg[1] = node->burgY[j];
				xdirseg[2] = node->burgZ[j];
				NormalizeVec(xdirseg);
				zdirseg[0] = node->nx[j];
				zdirseg[1] = node->ny[j];
				zdirseg[2] = node->nz[j];
				NormalizedCrossVector(zdirseg, xdirseg, ydirseg);
				Qseg[0][0] = xdirseg[0];
				Qseg[1][0] = xdirseg[1];
				Qseg[2][0] = xdirseg[2];
				Qseg[0][1] = ydirseg[0];
				Qseg[1][1] = ydirseg[1];
				Qseg[2][1] = ydirseg[2];
				Qseg[0][2] = zdirseg[0];
				Qseg[1][2] = zdirseg[1];
				Qseg[2][2] = zdirseg[2];

				Matrix33Transpose(Qseg, Qinvseg);

/*
 *				Transform everything into local coordinates for Jacobian
 *				calculation.
 */

				UnpackQMatrices(node, Q1, Qinv1);
				UnpackQMatrices(nbr, Q2, Qinv2);

				vectmp[0] = node->x;
				vectmp[1] = node->y;
				vectmp[2] = node->z;
				Matrix33Vector3Multiply(Qinvseg, vectmp, r1_loc);
				xtmp = nbr->x;
				ytmp = nbr->y;
				ztmp = nbr->z;
				PBCPOSITION(param, node->x, node->y, node->z, &xtmp, &ytmp, &ztmp);
				vectmp[0] = xtmp;
				vectmp[1] = ytmp;
				vectmp[2] = ztmp;
				Matrix33Vector3Multiply(Qinvseg, vectmp, r2_loc);

				dr[0] = r2_loc[0] - r1_loc[0];
				dr[1] = r2_loc[1] - r1_loc[1];
				dr[2] = r2_loc[2] - r1_loc[2];
				L = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				dr_norm[0] = dr[0]/L;
				dr_norm[1] = dr[1]/L;
				dr_norm[2] = dr[2]/L;

				vectmp[0] = node->burgX[j];
				vectmp[1] = node->burgY[j];
				vectmp[2] = node->burgZ[j];
				Matrix33Vector3Multiply(Qinvseg, vectmp, b_loc);

				Matrix33Mult33(Qinvseg, stressMat, tmpMat);
				Matrix33Mult33(tmpMat, Qseg, stressMatloc);
				
				vectmp[0] = node->fX;
				vectmp[1] = node->fY;
				vectmp[2] = node->fZ;
				Matrix33Vector3Multiply(Qinvseg, vectmp, f1_loc);

				vectmp[0] = nbr->fX;
				vectmp[1] = nbr->fY;
				vectmp[2] = nbr->fZ;
				Matrix33Vector3Multiply(Qinvseg, vectmp, f2_loc);

/*
 *				Calculate each of the partial Jacobians: core energy
 *				(line tension), length change, and applied stress. For
 *				each partial Jacobian, there are four sub-Jacobians 
 *				corresponding to 1-1, 1-2, 2-1, and 2-2 interactions.
 *				If either node is constrained to move in its local 
 *				x direction only (xonly# flag), then illegal terms need
 *				to be zeroed out using the ConstrainJacobian function.
 */

				B = 1/param->MobEdge;
				ndof1 = node->ndof;
				ndof2 = nbr->ndof;

/*
 *				Core energy partial Jacobian
 */
				if (param->elasticinteraction) {
					Ecore = param->Ecore;
				} else {
					Ecore = 0.5*param->TensionFactor*param->shearModulus;
				}

				bmag2 =	(node->burgX[j]*node->burgX[j] + node->burgY[j]*node->burgY[j] + 
								node->burgZ[j]*node->burgZ[j]);
				coreFactor = Ecore*bmag2/pow(L,5)*
						((1+param->pois)*dr[0]*dr[0]+(1-2*param->pois)*dr[1]*dr[1])/
						(1-param->pois);

				J0[0][0] = -coreFactor*dr[1]*dr[1];
				J0[0][1] = coreFactor*dr[0]*dr[1];
				J0[1][0] = coreFactor*dr[0]*dr[1];
				J0[1][1] = -coreFactor*dr[0]*dr[0];
				J0[0][2] = 0;
				J0[1][2] = 0;
				J0[2][0] = 0;
				J0[2][1] = 0;
				J0[2][2] = 0;

				locs[0] = 1; locs[1] = 1; locs[2] = 1; locs[3] = 1;

				ConstrainJacobian(Qseg, Qinvseg, Q1, Qinv1, Q2, Qinv2,
					J0, ndof1, ndof2, locs, J01a, J01b, J02a, J02b);
				
				for (i1=0; i1<3; i1++) {
					for (i2=0; i2<3; i2++) {
						Jcore[i1][i2] = J01a[i1][i2]/(B*Lsum1/2);
						Jcore[i1][i2+3] = -J01b[i1][i2]/(B*Lsum1/2);
						Jcore[i1+3][i2] = -J02a[i1][i2]/(B*Lsum2/2);
						Jcore[i1+3][i2+3] = J02b[i1][i2]/(B*Lsum2/2);
					}
				} 

/*
 *				Applied stress partial Jacobian
 */
				Matrix33Vector3Multiply(stressMatloc, b_loc, vectmp);
				cross(vectmp, dr_norm, fapp);
				fapp[2] = 0;
				cross(dr_norm, fapp, vectmp);
				s = vectmp[2];
				J0[0][0] = 0; J0[0][2] = 0; 
				J0[1][1] = 0; J0[1][2] = 0;
				J0[2][0] = 0; J0[2][1] = 0; J0[2][2] = 0;
				J0[0][1] = -0.5*s;
				J0[1][0] = 0.5*s;

				locs[0] = 0; locs[1] = 1; locs[2] = 0; locs[3] = 1;

				ConstrainJacobian(Qseg, Qinvseg, Q1, Qinv1, Q2, Qinv2,
					J0, ndof1, ndof2, locs, J01a, J01b, J02a, J02b);
				
				for (i1=0; i1<3; i1++) {
					for (i2=0; i2<3; i2++) {
						Japp[i1][i2] = J01a[i1][i2]/(B*Lsum1/2);
						Japp[i1][i2+3] = J01b[i1][i2]/(B*Lsum1/2);
						Japp[i1+3][i2] = -J02a[i1][i2]/(B*Lsum2/2);
						Japp[i1+3][i2+3] = J02b[i1][i2]/(B*Lsum2/2);
					}
				}

/*
 * 				Length change partial Jacobian
 */
				
				Matrix31Vector3Mult(f1_loc, dr_norm, tmpMat); 
				for (i1=0; i1<3; i1++) {
					for (i2=0; i2<3; i2++) {
						J0[i1][i2] = -tmpMat[i1][i2]/Lsum1;
					}
				}

				locs[0] = 1; locs[1] = 1; locs[2] = 0; locs[3] = 0;

				ConstrainJacobian(Qseg, Qinvseg, Q1, Qinv1, Q2, Qinv2,
					J0, ndof1, ndof2, locs, J01a, J01b, tmpMat, tmpMat2);
								
				
				Matrix31Vector3Mult(f2_loc, dr_norm, tmpMat); 
				for (i1=0; i1<3; i1++) {
					for (i2=0; i2<3; i2++) {
						J0[i1][i2] = tmpMat[i1][i2]/Lsum2;
					}
				}

				locs[0] = 0; locs[1] = 0; locs[2] = 1; locs[3] = 1;

				ConstrainJacobian(Qseg, Qinvseg, Q1, Qinv1, Q2, Qinv2,
					J0, ndof1, ndof2, locs, tmpMat, tmpMat2, J02a, J02b);
				
				for (i1=0; i1<3; i1++) {
					for (i2=0; i2<3; i2++) {
						JL[i1][i2] = -J01a[i1][i2]/(B*Lsum1/2);
						JL[i1][i2+3] = J01b[i1][i2]/(B*Lsum1/2);
						JL[i1+3][i2] = J02a[i1][i2]/(B*Lsum2/2);
						JL[i1+3][i2+3] = -J02b[i1][i2]/(B*Lsum2/2);
					}
				}
/*
 *				Add them all together to get the segment Jacobian.
 */
				for (i1=0; i1<6; i1++) {
					for (i2=0; i2<6; i2++) {
						Jseg[i1][i2] = Jcore[i1][i2] + Japp[i1][i2] + JL[i1][i2];
					}
				}

/*
 *				Add these contributions to the global Jacobian. If node2 isn't flagged,
 *				then we don't add its terms.
 */
				subJ[0][0] = Jseg[0][0]; subJ[0][1] = Jseg[0][1]; subJ[0][2] = Jseg[0][2];
				subJ[1][0] = Jseg[1][0]; subJ[1][1] = Jseg[1][1]; subJ[1][2] = Jseg[1][2];
				subJ[2][0] = Jseg[2][0]; subJ[2][1] = Jseg[2][1]; subJ[2][2] = Jseg[2][2];
				AddJterms(param, node, node, subJ, Jmob);

				subJ[0][0] = Jseg[0][3]; subJ[0][1] = Jseg[0][4]; subJ[0][2] = Jseg[0][5];
				subJ[1][0] = Jseg[1][3]; subJ[1][1] = Jseg[1][4]; subJ[1][2] = Jseg[1][5];
				subJ[2][0] = Jseg[2][3]; subJ[2][1] = Jseg[2][4]; subJ[2][2] = Jseg[2][5];
				AddJterms(param, node, nbr, subJ, Jmob);

				subJ[0][0] = Jseg[3][0]; subJ[0][1] = Jseg[3][1]; subJ[0][2] = Jseg[3][2];
				subJ[1][0] = Jseg[4][0]; subJ[1][1] = Jseg[4][1]; subJ[1][2] = Jseg[4][2];
				subJ[2][0] = Jseg[5][0]; subJ[2][1] = Jseg[5][1]; subJ[2][2] = Jseg[5][2];
				AddJterms(param, nbr, node, subJ, Jmob);

				if ((nbr->flags & UPDATE_NODE_J)!=0) {
					subJ[0][0] = Jseg[3][3]; subJ[0][1] = Jseg[3][4]; subJ[0][2] = Jseg[3][5];
					subJ[1][0] = Jseg[4][3]; subJ[1][1] = Jseg[4][4]; subJ[1][2] = Jseg[4][5];
					subJ[2][0] = Jseg[5][3]; subJ[2][1] = Jseg[5][4]; subJ[2][2] = Jseg[5][5];
					AddJterms(param, nbr, nbr, subJ, Jmob);
				}
			} //loop over all neighbors of ith node
		} //loop over all nodes

/*
 * 		Add all the elastic interaction terms.
 */
		for (eint=0; eint<esize; eint++) {
/*
 * 			Unpack nodal information.
 */
			node1 = enhancelist[eint][0];
			node2 = enhancelist[eint][1];
			node3 = enhancelist[eint][2];
			node4 = enhancelist[eint][3];
			int tag1 = node1->myTag.index;
			int tag2 = node2->myTag.index;
			int tag3 = node3->myTag.index;
			int tag4 = node4->myTag.index;
			Lsum1 = Lsumvec[node1->myTag.index];
			Lsum2 = Lsumvec[node2->myTag.index];
			Lsum3 = Lsumvec[node3->myTag.index];
			Lsum4 = Lsumvec[node4->myTag.index];

/*
 *			Calculate and add these terms to the Jacobian.
 */
			EnhanceJacobian(home, node1, node2, node3, node4,
							B, Lsum1, Lsum2, Lsum3, Lsum4, Jmob);
		}


/*
 * 		Unflag all nodes.
 */
		for (i=0; i < N; i++) {
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
            node->flags &= ~UPDATE_NODE_J;
		}

/*
 * 		Clean up memory.
 */
		free(Lsumvec);
		for(i=0;i<maxenhnum;i++){
			free(enhancelist[i]);
		}
		free(enhancelist);

		return;
}

/*------------------------------------------------------------------------
 *
 *      Function:    PreserveNodalData
 *      Description: Both old and new values for certain nodal
 *                   data items are required during timestep
 *                   integration and calculating plastic strain.
 *                   This function copies the values for specified
 *                   items into the appropriate variables.
 *                   
 *-----------------------------------------------------------------------*/
static void PreserveNodalData(Home_t *home)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
 
/*
 *          Save previous nodal position.
 */
			node->oldx = node->x;
			node->oldy = node->y;
			node->oldz = node->z;

/*
 *          Make a copy of the current velocity.
 */
			node->currvX = node->vX;    node->oldvX = node->vX;
			node->currvY = node->vY;    node->oldvY = node->vY;
			node->currvZ = node->vZ;    node->oldvZ = node->vZ;
        }

        node = home->ghostNodeQ;

        while (node != (Node_t *)NULL) {

			node->oldx = node->x;
			node->oldy = node->y;
			node->oldz = node->z;

			node->currvX = node->vX;    node->oldvX = node->vX;
			node->currvY = node->vY;    node->oldvY = node->vY;
			node->currvZ = node->vZ;    node->oldvZ = node->vZ;

            node = node->next;
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    ImplicitIntegrator
 *      Description: Implements a numerical timestep integrator using
 *                   the implicitly-solved Trapezoid integration method.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
static int Implicit_Step(Home_t *home , real8 newDT)
{
		int     i, j, convergent, maxIterations, tmp;
		int 	ind1, ind2, ind3;
		int     iter, globalIterError, mobIterError;
		int     dumpErrorData = 1, doAll = 1;
		int 	changedJ, errInc, count, NJ, Ndof;
		real8   errMax, globalErrMax, oldglobalErrMax, drn, relerrMax;
		real8   oldx, oldy, oldz, globalRelErrMax;
		real8   localVals[3], globalVals[3];
		real8 	Q[3][3], Qinv[3][3];
		real8	rtmp[3], rout[3], vtmp[3], vout[3];
		real8	thisx, thisy, thisz;
		real8 	errdX, errdY, errdZ, f1, f2, f3;
		real8	dx, dy, dz, xtmp, ytmp, ztmp;
		Node_t  *node;
		Param_t *param;

		param = home->param;
		int N = home->newNodeKeyPtr;

#ifndef _LIS
		Fatal("Linear Algebra Library (-D_LIS) is required for implicit integrator!");
                return (0);
#else
		LIS_VECTOR 	negf_lis, dr_lis;
		LIS_MATRIX 	J_lis, Jmob_lis;
		LIS_SOLVER 	solver;

/*
 * 		Allocate memory for arrays.
 */
		cs      *J, *Jmobtrp, *Jmob, *Eye, *tmp_cs;
		real8	*f, *r0, *v0, *Jval;
		int 	*Jrow, *Jcol;
		tmp = 3*(N+1);
		f   = (real8 *)malloc(tmp * sizeof(real8));
		r0  = (real8 *)malloc(tmp * sizeof(real8));
		v0  = (real8 *)malloc(tmp * sizeof(real8));

/*
 *		Set the initial guess to the current nodal positions
 */
		for (i=0; i < N; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if (node->ndof==0) continue;

			node->x = node->oldx;
			node->y = node->oldy;
			node->z = node->oldz;
		}
		
		mobIterError  = 0;
		errInc        = 0;
		globalVals[0] = 0.0;
		globalVals[1] = 0.0;
		globalVals[2] = 0.0;
		
/*
 *		Calculate nodal forces and velocities
 */
    	NodeForce(home, FULL);
        mobIterError = CalcNodeVelocities(home, 0, doAll);
        CommSendVelocity(home);
		if (mobIterError) return(0);

/*
 *      Loop until we converge on a time step.  
 */
        convergent    = 0;
        maxIterations = 20;

/*
 * 		Calculate coordinate systems and transformation matrices for every
 * 		node.
 */		
		for (i=0; i < N; i++) {
	        node = home->nodeKeys[i];
	        if (node == (Node_t *)NULL) continue;
			node->flags |= UPDATE_NODE_J;
		}
		Ndof = 0;
		AssignCoordSys(home, &Ndof, -1);
		
/*
 *		Calculate the Jacobian
 */ 
		Jmobtrp = cs_spalloc(Ndof, Ndof, 2*Ndof, 1, 1);
		UpdateJacobian(home,1,&changedJ,Jmobtrp);  // Jmobtrp gets updated here. 1 means update all nodes. changedJ is useless.
		
		Jmob = cs_compress(Jmobtrp);
		cs_spfree(Jmobtrp);
		cs_dupl(Jmob);

		tmp_cs = cs_spalloc(Ndof, Ndof, Ndof, 1, 1) ;    /* create triplet identity matrix */
		for (i = 0 ; i < Ndof ; i++) cs_entry(tmp_cs, i, i, 1) ;
		Eye = cs_compress(tmp_cs);
		cs_spfree(tmp_cs);

/*
 *		Calculate the iteration matrix. Diagonal entries first, then off-
 *		diagonal entries.
 */			
		lis_matrix_create(0,&J_lis);
		lis_matrix_set_size(J_lis,0,Ndof);

		J = cs_add(Eye, Jmob, 1, -newDT/2);

		int column = 0;
		for (j=0; j < J->nzmax; j++) {
			if (j == J->p[column+1]) column++;
			
			lis_matrix_set_value(0,J->i[j], column, J->x[j], J_lis);
		}
		lis_matrix_assemble(J_lis);

/*
 *		Set the initial guess to the current nodal positions.
 *		Convert the old position and velocity vectors to local coordinates
 */
		for (i=0; i < N; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			if ( node->ndof==0) continue;

			UnpackQMatrices(node, Q, Qinv);

			rtmp[0] = node->oldx;   vtmp[0] = node->currvX;
			rtmp[1] = node->oldy;   vtmp[1] = node->currvY;
			rtmp[2] = node->oldz;   vtmp[2] = node->currvZ;

			Matrix33Vector3Multiply(Qinv, rtmp, rout);
			Matrix33Vector3Multiply(Qinv, vtmp, vout);
			
			r0[i    ] = rout[0];   v0[i    ] = vout[0];
			r0[i+N  ] = rout[1];   v0[i+N  ] = vout[1];
			r0[i+2*N] = rout[2];   v0[i+2*N] = vout[2];
		}
		
		for (iter = 0; iter < maxIterations; iter++) {
/*
 *			Recalculate nodal force and velocity at the new positions
 */
			if (iter>0) { //No need to update velocities on 1st iteration
				NodeForce(home, FULL);
				mobIterError = CalcNodeVelocities(home, 0, doAll);
				CommSendVelocity(home);
				if (mobIterError) return(0);
			}

/*
 *			Calculate the residual vector and maximum error in local
 *			coordinates.
 */
			errMax    = 0.0;
			relerrMax = 0.0;

			for ( i=0 ; i < N ; i++) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				if ( node->ndof==0) continue;

				UnpackQMatrices(node, Q, Qinv);

				thisx = node->x;
				thisy = node->y;
				thisz = node->z;
				PBCPOSITION(param, node->oldx, node->oldy, node->oldz,
				                   &thisx    , &thisy    , &thisz   );

				rtmp[0] = thisx;
				rtmp[1] = thisy;
				rtmp[2] = thisz;					
				vtmp[0] = node->vX;
				vtmp[1] = node->vY;
				vtmp[2] = node->vZ;

				Matrix33Vector3Multiply(Qinv, rtmp, rout);
				Matrix33Vector3Multiply(Qinv, vtmp, vout);

				if (node->dofx != 0) {
					f1 = rout[0] - r0[i] - newDT/2*(vout[0] + v0[i]);
					f[node->dofx-1] = f1;
				} else {
					f1 = 0;
				}
				if (node->dofy != 0) {
					f2 = rout[1] - r0[i+N] - newDT/2*(vout[1] + v0[i+N]);
					f[node->dofy-1] = f2;
				} else {
					f2 = 0;
				}
				if (node->dofz != 0) {
					f3 = rout[2] - r0[i+2*N] - newDT/2*(vout[2] + v0[i+2*N]);
					f[node->dofz-1] = f3;
				} else {
					f3 = 0;
				}
				
				real8 thiserr = sqrt( f1*f1 + f2*f2 + f3*f3);
				errMax = MAX(errMax, thiserr);
				
				//Calculate relative error
				errdX = thisx - node->oldx;
				errdY = thisy - node->oldy;
				errdZ = thisz - node->oldz;
				drn   = sqrt( errdX*errdX + errdY*errdY + errdZ*errdZ );
				
				if (thiserr > param->rTolth) {
					if (drn > param->rTolth / param->rTolrel) {
						relerrMax = MAX( relerrMax , thiserr/drn);
					} else {
						relerrMax = 2.0*param->rTolrel;
					}
				}
				
			/*	if ((thiserr < param->rTol) && (thiserr < param->rTolth || thiserr/drn < param->rTolrel) && iter>0) {
					node->subgroup = 0;
				} else {
					node->subgroup = 1;
				}*/
			}

/*
 *			Need to find largest errMax from among all domains.
 */
#if PARALLEL
			localVals[0] = errMax;
			localVals[1] = (real8)mobIterError;
			localVals[2] = relerrMax;

			MPI_Allreduce(localVals, globalVals, 3, MPI_DOUBLE, MPI_MAX, 
			              MPI_COMM_WORLD);

			globalErrMax    = globalVals[0];
			globalIterError = globalVals[1];
			globalRelErrMax = globalVals[2];
#else
			globalErrMax    = errMax;
			globalIterError = mobIterError;
			globalRelErrMax = relerrMax;
#endif

/*
 *			If the error is within the tolerance, we've reached
 *			convergence so we can accept this deltaT.  Otherwise
 *			continue iterating, unless max number of iterations has reached.
 *			If it has, exit the for loop and cut the time step.  
 *			Note: we need to reposition both local nodes and ghost nodes!
 */
			if (globalErrMax < param->rTol && globalRelErrMax < param->rTolrel
			                               && iter > 0) {
				convergent = 1;
				break;
			} else {
				if (iter == maxIterations-1) continue;
				
				if (iter>0) {
					if (globalErrMax > oldglobalErrMax) {
						errInc++;
						if (errInc>=3) break;
					}
				}
				oldglobalErrMax = globalErrMax;

/*				Solve for the new nodal positions and velocities using the 
 *				Newton-Raphson method.
 */
				//Instantiate lis vectors and solver.
				lis_vector_create(0, &negf_lis);
				lis_vector_create(0, &dr_lis);

				lis_solver_create(&solver);
				lis_solver_set_option("-i bicg -p none",solver);
				lis_solver_set_option("-tol 1.0e-12",solver);

				lis_vector_set_size(negf_lis, Ndof, 0);
				lis_vector_set_size(dr_lis, Ndof, 0);

				//Set values of negf_lis vector.
				for (i=0; i < N; i++) {
					if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
					if ( node->ndof==0) continue;

					ind1 = node->dofx;
					ind2 = node->dofy;
					ind3 = node->dofz;

					if (ind1 != 0) lis_vector_set_value(0, ind1-1, -f[ind1-1], negf_lis);
					if (ind2 != 0) lis_vector_set_value(0, ind2-1, -f[ind2-1], negf_lis);
					if (ind3 != 0) lis_vector_set_value(0, ind3-1, -f[ind3-1], negf_lis);
				}

				//Solve the system
				lis_solve(J_lis,negf_lis,dr_lis,solver);

				//Extract the new nodal positions and rotate to global coords
				for (i=0; i < N; i++) {
					if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
					if ( node->ndof==0) continue;

					ind1 = node->dofx;
					ind2 = node->dofy;
					ind3 = node->dofz;

					if (ind1 != 0) lis_vector_get_value(dr_lis, ind1-1, &dx);
					else           dx = 0;
					
					if (ind2 != 0) lis_vector_get_value(dr_lis, ind2-1, &dy);
					else           dy = 0;
					
					if (ind3 != 0) lis_vector_get_value(dr_lis, ind3-1, &dz);
					else           dz = 0;

					UnpackQMatrices(node, Q, Qinv);

					rtmp[0] = dx;
					rtmp[1] = dy;
					rtmp[2] = dz;					

					Matrix33Vector3Multiply(Q, rtmp, rout);

					xtmp = node->x + rout[0];
					ytmp = node->y + rout[1];
					ztmp = node->z + rout[2];

					FoldBox(param,&xtmp,&ytmp,&ztmp);

					node->x = xtmp;
					node->y = ytmp;
					node->z = ztmp;					
				}

				//Kill lis vectors.
				lis_solver_destroy(solver);
				lis_vector_destroy(negf_lis);
				lis_vector_destroy(dr_lis);

			}  /* not convergent */
		}  /* for (iter = 0; ...) */

/*
 *		Need to kill the lis object because we will reinstantiate it.
 */
		lis_matrix_destroy(J_lis);

/*
 *		Clean up memory.
 */
		free(r0);
		free(v0);
		free(f);
		cs_spfree(Jmob);
		cs_spfree(J) ;
		cs_spfree(Eye); 

		return(convergent);
#endif
}


/*------------------------------------------------------------------------
 *
 *      Function:     cutTimeStepSize
 *      Description:  This function cuts the delta T by a configured 
 *                    factor, if we have not reached convergence
 *
 *-----------------------------------------------------------------------*/
static void cutTimeStepSize(Home_t *home , real8 *newDT , int *incrDelta)
{
		Param_t *param;
		param = home->param;

		*newDT         *=  param->dtDecrementFact;
		*incrDelta      =  0;
		param->deltaTT  = *newDT;

		if ((*newDT < 1.0e-20) && (home->myDomain == 0)) {
			Fatal("TrapezoidIntegrator(): Timestep has dropped below\n"
				  "minimal threshold to %e.  Aborting!", *newDT);
		}
}


/*------------------------------------------------------------------------
 *
 *      Function:    ImplicitIntegrator
 *      Description: Implements a numerical timestep integrator using
 *                   the implicitly-solved Trapezoid integration method.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
void ImplicitIntegrator(Home_t *home)
{
		int     i, convergent, incrDelta;
		int     globalIterError;
		int     doAll = 1, N = home->newNodeKeyPtr;
		real8   errMax, localVals[2], globalVals[2];
		real8   newDT , oldx  , oldy  , oldz     , drn;
		real8	errx  , erry  , errz  , thiserr  , globalErrMax;
		real8   errdX , errdY , errdZ , relerrMax, globalRelErrMax;
		real8   tmp1  , tmp2  , tmp3  , tmp4     , factor;
		Node_t  *node;
		Param_t *param;

		param = home->param;
		
		int    *col = malloc(N * sizeof(int     ));
		real8  **X0 = malloc(N * sizeof(real8 * ));
		real8  **X1 = malloc(N * sizeof(real8 * ));
		real8  **X2 = malloc(N * sizeof(real8 * ));
		
		for ( i = 0 ; i < N ; i++ ) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			X0[i] = malloc( 3 * sizeof(real8 ));
			X1[i] = malloc( 3 * sizeof(real8 ));
			X2[i] = malloc( 3 * sizeof(real8 ));
			
			X0[i][0] = node->x;
			X0[i][1] = node->y;
			X0[i][2] = node->z;
		}

/*
 * 		Set the initial time step.
 */
		newDT = MIN(param->maxDT, param->nextDT);
		if (newDT <= 0.0) newDT = param->maxDT;
		param->deltaTT = newDT;

/*
 *      Preserve certain items of nodal data for later use.
 */
        PreserveNodalData(home);

/*
 *      Loop until we converge on a time step.  First step is to
 *      use the current positions and velocities and reposition the
 *      nodes to where they would be after the suggested delta time.
 */
        convergent = 0;
        incrDelta  = 1;

        int iTry1=0;
		int iTry2=0;
		
		while (!convergent) {
			errMax        = 0.0;
			relerrMax     = 0.0;
			globalVals[0] = 0.0;
			globalVals[1] = 0.0;
		
/*
 *			Taking one step of size dt and storing the positions in X1
 */
			convergent = Implicit_Step(home , newDT);
			if (!convergent) {
				iTry1++;
			/*	if (iTry1==1) {
					int ioGroup = home->ioGroupNum;
					char debug_data [50];
					sprintf (debug_data, "first%d", home->cycle);
					WriteAtomEye(home, debug_data, ioGroup, 1, 1, 1);
				}*/
			
				cutTimeStepSize(home , &newDT , &incrDelta);
				continue;
			}
		//	continue;
			
			iTry2++;
			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				X1[i][0] = node->x;
				X1[i][1] = node->y;
				X1[i][2] = node->z;
			}
			
/*
 *			Taking two RK steps of size dt/2
 *			First take one step of size dt/2 starting form the old positions
 */
			convergent = Implicit_Step(home , newDT*0.5);
			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				node->oldx = node->x;   node->currvX = node->vX;
				node->oldy = node->y;   node->currvY = node->vY;
				node->oldz = node->z;   node->currvZ = node->vZ;
			}
			
/*
 *			Now take the second half step and evaluate the truncation error
 */
			convergent = Implicit_Step(home , newDT*0.5);
			
			for ( i = 0 ; i < N ; i++ ) {
				if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
				
				node->oldx = X0[i][0];  X2[i][0] = node->x;
				node->oldy = X0[i][1];  X2[i][1] = node->y;
				node->oldz = X0[i][2];  X2[i][2] = node->z;
				
				node->currvX = node->oldvX;
				node->currvY = node->oldvY;
				node->currvZ = node->oldvZ;
				
/*
 *				Calculate the truncation error on this node
 */
				errx = X2[i][0] - X1[i][0];
				erry = X2[i][1] - X1[i][1];
				errz = X2[i][2] - X1[i][2];
				ZImage(param, &errx, &erry, &errz);
				
				thiserr = sqrt( errx*errx + erry*erry + errz*errz );
				errMax  = MAX ( errMax , thiserr);

/*
 *				Calculate the relative error on this node
 */
				errdX = node->x - X0[i][0];
				errdY = node->y - X0[i][1];
				errdZ = node->z - X0[i][2];
				ZImage(param, &errdX, &errdY, &errdZ);
				drn   = sqrt( errdX*errdX + errdY*errdY + errdZ*errdZ );

				if (thiserr > param->rTolth ) {
					if (drn > param->rTolth / param->rTolrel ) {
						relerrMax = MAX( relerrMax , thiserr/drn );
					} else {
						relerrMax = 2.0*param->rTolrel;
					}
				}
				
				if (iTry2==1) {
					if ((thiserr < param->rTol) && (thiserr < param->rTolth || thiserr/drn < param->rTolrel)) {
						col[i] = 0;
					} else {
						col[i] = 1;
					}
				}
			}
			
/*
 *			Need to find largest errMax from among all domains.
 */
#if PARALLEL
			localVals[0] = errMax;
			localVals[1] = relerrMax;

			MPI_Allreduce(localVals, globalVals, 2, MPI_DOUBLE, MPI_MAX, 
			              MPI_COMM_WORLD);

			globalErrMax    = globalVals[0];
			globalRelErrMax = globalVals[1];
#else
			globalErrMax    = errMax;
			globalRelErrMax = relerrMax;
#endif

/*
 *			If the error is within the tolerance, we've reached
 *			convergence so we can accept this deltaT.  Otherwise
 *			reposition the nodes and try again.  Note: we need to
 *			reposition both local nodes and ghost nodes!
 */
			if (globalErrMax < param->rTol && globalRelErrMax < param->rTolrel) {
				convergent = 1;
				break;
			} else {
				convergent = 0;
				cutTimeStepSize(home , &newDT , &incrDelta);
			}
		}
		
	/*	for (i=0 ; i<N ; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			node->subgroup = col[i];
		}*/
		
	/*	int ioGroup = home->ioGroupNum;
		char debug_data [50];
		sprintf (debug_data, "second%d", home->cycle);
		WriteAtomEye(home, debug_data, ioGroup, 1, 1, 1);*/

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
        param->deltaTT   = newDT;
        param->realdt    = newDT;
        param->timeStart = param->timeNow;

        if (incrDelta) {
            if (param->dtVariableAdjustment) {
                tmp1 = pow(param->dtIncrementFact, param->dtExponent);
                tmp2 = globalErrMax/param->rTol;
                tmp3 = 1.0 / param->dtExponent;
                tmp4 = pow(1.0/(1.0+(tmp1-1.0)*tmp2), tmp3);
                factor = param->dtIncrementFact * tmp4;
				
                param->nextDT = MIN(param->maxDT, newDT*factor);
            } else {
                param->nextDT = MIN(param->maxDT, newDT*param->dtIncrementFact);
            }
        } else {
            param->nextDT = newDT;
        }
	//printf("iTry= %d  %d\n" , iTry1 , iTry2);
		
/*
 *		Freeing memory
 */ 
		for (i=0; i < home->newNodeKeyPtr; i++) {
			if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
			
			free( X0[i] );
			free( X1[i] );
			free( X2[i] );
		}
		free( X0 );
		free( X1 );
		free( X2 );
		free( col);
		
		return;
}

#endif //_IMPLICIT
