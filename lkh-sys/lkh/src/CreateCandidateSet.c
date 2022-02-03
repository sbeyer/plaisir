#include "LKH.h"

/*
 * The CreateCandidateSet function determines for each node its set of incident
 * candidate edges.
 *
 * The Ascent function is called to determine a lower bound on the optimal tour 
 * using subgradient optimization. But only if the penalties (the Pi-values) is
 * not available on file. In the latter case, the penalties is read from the 
 * file, and the lower bound is computed from a minimum 1-tree.      
 *
 * The function GenerateCandidates is called to compute the Alpha-values and to 
 * associate to each node a set of incident candidate edges.  
 *
 * The CreateCandidateSet function itself is called from LKHmain.
 */

void CreateCandidateSet()
{
    GainType Cost, MaxAlpha;
    Node *Na;
    double EntryTime = GetTime();

    Norm = 9999;
    if ((MaxTrials == 0 &&
             (FirstNode->InitialSuc || InitialTourAlgorithm == SIERPINSKI ||
              InitialTourAlgorithm == MOORE))) {
        AddTourCandidates();
        goto End_CreateCandidateSet;
    }
    if (TraceLevel >= 2)
        printff("Creating candidates ...\n");
    if (MaxCandidates > 0 &&
            (CandidateSetType == QUADRANT || CandidateSetType == NN)) {
        AddTourCandidates();
        if (CandidateSetSymmetric)
            SymmetrizeCandidateSet();
        goto End_CreateCandidateSet;
    }

    Na = FirstNode;
    do
        Na->Pi = 0;
    while ((Na = Na->Suc) != FirstNode);
    Cost = Ascent();

    LowerBound = (double) Cost / Precision;
    if (TraceLevel >= 1) {
        printff("Lower bound = %0.1f", LowerBound);
        printff(", Ascent time = %0.2f sec.\n",
                fabs(GetTime() - EntryTime));
    }
    MaxAlpha = (GainType) fabs(Excess * Cost);
    if (CandidateSetType == DELAUNAY ||
            MaxCandidates == 0)
        OrderCandidateSet(MaxCandidates, MaxAlpha, CandidateSetSymmetric);
    else
        GenerateCandidates(MaxCandidates, MaxAlpha, CandidateSetSymmetric);

End_CreateCandidateSet:
    if (ExtraCandidates > 0) {
        AddExtraCandidates(ExtraCandidates,
                ExtraCandidateSetType,
                ExtraCandidateSetSymmetric);
        AddTourCandidates();
    }
    ResetCandidateSet();
    if (MaxTrials > 0 ||
            (InitialTourAlgorithm != SIERPINSKI &&
             InitialTourAlgorithm != MOORE)) {
        Na = FirstNode;
        do {
            if (!Na->CandidateSet || !Na->CandidateSet[0].To) {
                if (MaxCandidates == 0)
                    eprintf
                        ("MAX_CANDIDATES = 0: Node %d has no candidates",
                         Na->Id);
                else
                    eprintf("Node %d has no candidates", Na->Id);
            }
        }
        while ((Na = Na->Suc) != FirstNode);
    }
    if (TraceLevel >= 1) {
        CandidateReport();
        printff("Preprocessing time = %0.2f sec.\n",
                fabs(GetTime() - EntryTime));
    }
}
