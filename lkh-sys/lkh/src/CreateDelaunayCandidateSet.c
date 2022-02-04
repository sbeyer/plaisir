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
        if (TraceLevel >= 2)
            printff("done\n");
    }
}
