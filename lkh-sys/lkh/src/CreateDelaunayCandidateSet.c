#include "LKH.h"
#include "Delaunay.h"

/*
 * The CreateDelaunayCandidateSet function determines for each node its set
 * of incident candidate edges. The edges are found by Delaunay triangulation.
 *
 * The function is called from CreateCandidateSet.
 */

static int Level = 0;

void CreateDelaunayCandidateSet()
{
    Node *From, *To;
    point *u, *v;
    edge *e_start, *e;
    int d, i, Count;

    if (TraceLevel >= 2)
        printff("Creating Delaunay candidate set ... ");
    if (Level == 0 && MaxCandidates == 0) {
        AddTourCandidates();
        From = FirstNode;
        do {
            if (!From->CandidateSet)
                eprintf("MAX_CANDIDATES = 0: No candidates");
        } while ((From = From->Suc) != FirstNode);
        if (TraceLevel >= 2)
            printff("done\n");
        return;
    }

    /* Find the Delaunay edges */
    delaunay(Dimension);

    /* Add the Delaunay edges to the candidate set */
    for (i = 0; i < Dimension; i++) {
        u = &p_array[i];
        From = &NodeSet[u->id];
        e_start = e = u->entry_pt;
        Count = 0;
        do {
            v = Other_point(e, u);
            if (u < v) {
                To = &NodeSet[v->id];
                d = D(From, To);
                AddCandidate(From, To, d, 1);
                AddCandidate(To, From, d, 1);
            }
        } while ((e = Next(e, u)) != e_start && ++Count < Dimension);
    }
    free_memory();
    if (Level == 0) {
        AddTourCandidates();
        if (ExtraCandidates < 2) {
            /* Add quadrant neighbors if any node has less than two candidates.
               That is, if it should happen that delaunay_edges fails. */
            From = FirstNode;
            do {
                if (From->CandidateSet == 0 ||
                    From->CandidateSet[0].To == 0 ||
                    From->CandidateSet[1].To == 0) {
                    if (TraceLevel >= 2)
                        printff("*** Not complete ***\n");
                    AddExtraCandidates(4,
                                       QUADRANT, 1);
                    break;
                }
            } while ((From = From->Suc) != FirstNode);
        }
        if (TraceLevel >= 2)
            printff("done\n");
    }
}
