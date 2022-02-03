#include "LKH.h"

/*
 * The FixedOrCommonCandidates function returns the number of fixed or
 * common candidate edges emanating from a given node, N.
 */

int FixedOrCommonCandidates(Node * N)
{
    int Count = 0;

    Count = N->FixedTo2 ? 2 : N->FixedTo1 ? 1 : 0;
    if (Count > 2)
        eprintf("Node %d has more than two required candidate edges",
                N->Id);
    return Count;
}
