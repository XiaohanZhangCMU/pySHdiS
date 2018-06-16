/**************************************************************************
 *
 *      Module:      CommSendVelocitySub.c
 *      Description: Contains functions necessary for communicating
 *                   updated nodal velocity data between domains
 *                   that own the nodes and domains that have the nodes
 *                   exported as ghost nodes.
 *
 *                   This function is intended for use with the 
 *                   subcycle integrator and only nodes with their
 *                   CommSend flag toggled will be sent to other domains.
 *
 *      Includes functions:
 *
 *          CommSendVelocitySub()
 *          PackVelocity()
 *          Unpackvelocity()
 *
 ***************************************************************************/

#include "Home.h"
#include "Comm.h"

#define VEL_FLTS_PER_NODE 4
#define VEL_FLTS_PER_ARM  5
#define VEL_FLTS_EXTRA    1


/*------------------------------------------------------------------------
 *
 *      Function:    PackVelocity
 *      Description: Allocate and pack data containing
 *                   velocity data for local nodes that are exported
 *                   as ghost nodes to neighboring domains.
 *
 *-----------------------------------------------------------------------*/
static void PackVelocity(Home_t *home, int subGroup)
{
        int            i, j, k, bufIndex, sgroup;
        int            domainIdx, cellIndex, nodesSent, armsSent, valCount;
        real8          *outBuf;
        Cell_t         *cell;
        Node_t         *node;
        RemoteDomain_t *remDom;
		
		sgroup = 0;
		if (subGroup >= GROUP0) sgroup = subGroup - 10;
/*
 *      Loop over neighboring domain.  Pack into the buffer for
 *      the domain, data for each local node that is contained
 *      in any cell exported from this domain to the neighboring
 *      domain.
 */
        for (i = 0; i < home->remoteDomainCount; i++) {

            domainIdx = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[domainIdx];

            nodesSent = 0;
            armsSent  = 0;

/*
 *          In order to know how large a buffer to allocate, we need
 *          to do a quick loop to total up the number of nodes and
 *          arms for which data will be sent to this neighboring
 *          domain.
 */
			for (j = 0; j < remDom->numExpCells; j++) {

				cellIndex = remDom->expCells[j];
				cell = home->cellKeys[cellIndex];

				node = cell->nodeQ;
				while (node != (Node_t *)NULL) {
					if ((node->myTag.domainID == home->myDomain)
                        && (node->CommSend[sgroup] == 1 || subGroup == FULL)) {
						nodesSent++;
					}
					node = node->nextInCell;
				}
			}

/*
 *          Allocate the buffer to be sent and pack it with data
 */
			valCount = (nodesSent * VEL_FLTS_PER_NODE) + VEL_FLTS_EXTRA;

            remDom->outBufLen = valCount * sizeof(real8);
            outBuf = (real8 *)malloc(remDom->outBufLen);
            bufIndex = 0;
            outBuf[bufIndex++] = nodesSent;

            for (j = 0; j < remDom->numExpCells; j++) {

                cellIndex = remDom->expCells[j];
                cell = home->cellKeys[cellIndex];
           
                node = cell->nodeQ;

                while (node != (Node_t *)NULL) {

                    if ( (node->myTag.domainID != home->myDomain) 
                         || (node->CommSend[sgroup] == 0 && subGroup != FULL)) {
						node = node->nextInCell;
                        continue;
                    }

                    outBuf[bufIndex++] = (real8)node->myTag.index;
					outBuf[bufIndex++] = node->vX;
					outBuf[bufIndex++] = node->vY;
					outBuf[bufIndex++] = node->vZ;

                    node = node->nextInCell;

                }  /* loop over nodes in cell */
            }  /* loop over exported cells */

            remDom->outBuf = (char *)outBuf;

/*
 *          We probably don't need this, but keep in in as a
 *          quick sanity check for now.
 */
            if (bufIndex > valCount) {
                Fatal("PackVelocity: packed %d values into %d element buffer",
                      bufIndex, valCount);
            }

        }  /* loop over remote domain */

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    UnpackVelocity
 *      Description: Copy velocity data for ghost nodes out
 *                   of the buffers sent from neighboring domains.
 *
 *-----------------------------------------------------------------------*/
static void UnpackVelocity(Home_t *home)
{
#ifdef PARALLEL
        int            i, j, k, arm, domIndex, bufIndex, numNodes, numArms;
        real8          fx, fy, fz;
        real8          *inBuf;
        Tag_t          tag, nbrTag;
        Node_t         *node;
        RemoteDomain_t *remDom;

        for (i = 0; i < home->remoteDomainCount; i++) {

            domIndex = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[domIndex];
            inBuf = (real8 *)remDom->inBuf;

            tag.domainID = domIndex;
            bufIndex = 0;
            numNodes = (int)inBuf[bufIndex++];

            for (j = 0; j < numNodes; j++) {

                tag.index = (int)inBuf[bufIndex++];

                if ((node = GetNodeFromTag(home, tag)) == (Node_t *)NULL) {
                    Fatal("UnpackVelocity: Remote node (%d,%d) is not a ghost",
                          tag.domainID, tag.index);
                }

                node->vX = inBuf[bufIndex++];
                node->vY = inBuf[bufIndex++];
                node->vZ = inBuf[bufIndex++];

            }  /* loop over all nodes in the buffer */
        }
#endif

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    CommSendVelocitySub
 *      Description: Driver function to send velocity data
 *                   for local nodes to neighboring domains and
 *                   receiving similar data for remote nodes this
 *                   domain maintains as ghost nodes.
 *
 *                   Intended for use with the subcycle integrator.
 *
 *-----------------------------------------------------------------------*/
void CommSendVelocitySub(Home_t *home, int subGroup)
{
#ifdef PARALLEL
        int            i, domainIdx;
        int            valCount, localBuffers = 0;
        RemoteDomain_t *remDom;
		
		if (home->numDomains == 1) return;

        TimerStart(home, COMM_SEND_VELOCITYSUB);
/*
 *      Pre-issue receives of message lengths from each neighbor
 */
        for (i = 0; i < home->remoteDomainCount; i++) {

            domainIdx = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Irecv(&remDom->inBufLen, 1, MPI_INT, domainIdx,
                      MSG_VELOCITY_LEN, MPI_COMM_WORLD, &home->inRequests[i]);
        }

/*
 *      Package up velocity data for the neighboring domains and send
 *      out the buffer lengths
 */
        PackVelocity(home, subGroup);

        for (i = 0; i < home->remoteDomainCount; i++) {

            domainIdx = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Isend(&remDom->outBufLen, 1, MPI_INT, domainIdx,
                      MSG_VELOCITY_LEN, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }


/*
 *      Wait for the length sends/receives to complete
 */
        MPI_Waitall(home->remoteDomainCount,home->outRequests, home->outStatus);
        MPI_Waitall(home->remoteDomainCount,home->inRequests, home->inStatus);
    
/*
 *      Allocate appropriately sized buffers for the incoming messages
 *      and pre-issue receives for buffers from all neighborning domains
 */
        for (i = 0; i < home->remoteDomainCount; i++) {

            domainIdx = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[domainIdx];

            valCount = remDom->inBufLen / sizeof(real8);
            remDom->inBuf = (char *)malloc(remDom->inBufLen);
            MPI_Irecv(remDom->inBuf, valCount, MPI_DOUBLE, domainIdx,
                      MSG_VELOCITY, MPI_COMM_WORLD, &home->inRequests[i]);
        }

/*
 *      Send local velocities out to all the neighboring domains.
 */
        for (i = 0; i < home->remoteDomainCount; i++) {

            domainIdx = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[domainIdx];

            valCount = remDom->outBufLen / sizeof(real8);
            MPI_Isend(remDom->outBuf, valCount, MPI_DOUBLE,
                      domainIdx, MSG_VELOCITY, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }

/*
 *      Wait for the data buffer sends/receives to complete and unpack
 *      the data.
 */
        MPI_Waitall(home->remoteDomainCount,home->outRequests, home->outStatus);
        MPI_Waitall(home->remoteDomainCount,home->inRequests, home->inStatus);

        UnpackVelocity(home);

/*
 *      Release all the message buffers...
 */
        for (i = 0; i < home->remoteDomainCount; i++) {

            domainIdx = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[domainIdx];

            localBuffers += remDom->inBufLen;
            localBuffers += remDom->outBufLen;

            free(remDom->inBuf);
            free(remDom->outBuf);
        }

/*
 *      Just some debug code to print the maximum buffer space
 *      used by any domain during the node velocity communications
 */
#if 0
{
        int globalBuffers = 0;

        MPI_Allreduce(&localBuffers, &globalBuffers, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);

        if (globalBuffers == localBuffers) {
            printf("  Task %d: Velocity comm total buffers = %dKb\n",
                   home->myDomain, globalBuffers / 1000);
        }
}
#endif

        TimerStop(home, COMM_SEND_VELOCITYSUB);

#endif /* if PARALLEL */

        return;
}
