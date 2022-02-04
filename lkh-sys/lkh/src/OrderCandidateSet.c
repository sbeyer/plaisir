#include "LKH.h"

/*
 * The OrderCandidateSet function augments the candidate set by using
 * transitive relatations in the following way. If the edges (i,j) and (j,k)
 * are contained the candidate set, then the edge (i,k) is added to the
 * candidate set. The alpha-value of each candidate edge is computed, and
 * the candidate edges associated with each node are ordered according to
 * their Alpha-values.
 *
 * The parameter MaxCandidates specifies the maximum number of candidate
 * edges allowed for each node.
 *
 * A non-zero value of Symmetric specifies that the candidate set is to be
 * complemented such that every candidate edge is associated with both its
 * two end nodes (in this way MaxCandidates may be exceeded).
 */

#define Ancestor OldPred        /* Nearest possible least ancestor */
#define AncestorSon OldSuc      /* Nearest son of Ancestor         */
#define OriginalCandidates Sons /* Number of original candidates   */
#define AlphaComputed OldSucExcluded    /* Have Alpha-values been computed? */
#define Level V

static int BetaValue(Node * From, Node * To);
static Candidate *FindCandidate(Node * From, Node * To);
#undef max
static int max(const int a, const int b);

void OrderCandidateSet(int MaxCandidates, int Symmetric)
{
    Node *From, *To;
    Candidate *NFrom, *NN;
    int Beta;

    if (TraceLevel >= 2)
        printff("Ordering candidates ... ");
    /* Add edges from the 1-tree to the candidate set */
    if (MaxCandidates > 0) {
        From = FirstNode;
        do {
            if ((To = From->Dad)) {
                AddCandidate(From, To, From->Cost, 0);
                AddCandidate(To, From, From->Cost, 0);
            }
        }
        while ((From = From->Suc) != FirstNode);
        AddCandidate(FirstNode, FirstNode->Next, FirstNode->NextCost, 0);
        AddCandidate(FirstNode->Next, FirstNode, FirstNode->NextCost, 0);
    }

    From = FirstNode;
    do {
        From->AlphaComputed = 0;
        From->Sons = 0;
    } while ((From = From->Suc) != FirstNode);
    From = FirstNode->Suc;
    From->Level = 0;
    From->Ancestor = From->AncestorSon = From;
    From->Beta = INT_MIN;
    From = From->Suc;
    do {
        From->Ancestor = To = From->Dad;
        From->Beta = !FixedOrCommon(From, To) ? From->Cost : INT_MIN;
        From->Level = To->Level + 1;
        To->Sons++;
        From->AncestorSon = From;
    }
    while ((From = From->Suc) != FirstNode);

    From = FirstNode->Suc->Suc;
    do {
        To = From->Dad;
        if (To->Sons == 1) {
            From->Beta = max(To->Beta, From->Beta);
            From->Ancestor = To->Ancestor;
            From->AncestorSon = To->AncestorSon;
        }
    }
    while ((From = From->Suc) != FirstNode);

    /* Compute Alpha-values for candidates */
    do {
        for (NFrom = From->CandidateSet; NFrom && (To = NFrom->To);
             NFrom++) {
            if (FixedOrCommon(From, To))
                NFrom->Alpha = INT_MIN;
            else if (From->FixedTo2 || To->FixedTo2)
                NFrom->Alpha = INT_MAX;
            else if (To->AlphaComputed && (NN = FindCandidate(To, From)))
                NFrom->Alpha = NN->Alpha;
            else {
                Beta = BetaValue(From, To);
                NFrom->Alpha =
                    Beta != INT_MIN ? max(NFrom->Cost - Beta, 0) : INT_MAX;
            }
        }
        From->AlphaComputed = 1;
    }
    while ((From = From->Suc) != FirstNode);

    /* Order candidates according to their Alpha-values */
    ResetCandidateSet();
    if (MaxCandidates > 0)
        TrimCandidateSet(MaxCandidates);
    AddTourCandidates();
    if (Symmetric)
        SymmetrizeCandidateSet();
    if (TraceLevel >= 2)
        printff("done\n");
}

/*
 * The BetaValue function computes the largest edge cost on the path
 * between two given nodes in the minimum spanning tree.
 */

static int BetaValue(Node * From, Node * To)
{
    Node *N1 = From, *N2 = To;
    int Beta = INT_MIN;

    if (To == From->Dad)
        return From->Cost;
    if (From == To->Dad)
        return To->Cost;
    if (From == FirstNode || To == FirstNode)
        return FirstNode->NextCost;

    /* Go upwards in the tree until the least common ancestor is met */
    while (N1->Ancestor != N2->Ancestor) {
        if (N1->Level > N2->Level) {
            if (N1->Beta > Beta)
                Beta = N1->Beta;
            N1 = N1->Ancestor;
        } else {
            if (N2->Beta > Beta)
                Beta = N2->Beta;
            N2 = N2->Ancestor;
        }
    }
    if (N1 == N2)
        return Beta;
    if (N1->AncestorSon != N2->AncestorSon)
        return max(Beta, max(N1->Beta, N2->Beta));
    if (N1->Level < N2->Level) {
        Node *t = N1;
        N1 = N2;
        N2 = t;
    }
    if (N1->Beta > N2->Beta)
        return max(Beta, N1->Beta);
    while (N1 != N2) {
        if (N1->Cost > Beta)
            Beta = N1->Cost;
        N1 = N1->Dad;
    }
    return Beta;
}

/*
 * The FindCandidate function returns the Candidate structure that is
 * associated with the node From and is pointing to the node To. The
 * function returns 0 if the search fails.
 */

static Candidate *FindCandidate(Node * From, Node * To)
{
    Candidate *NFrom;
    for (NFrom = From->CandidateSet; NFrom->To; NFrom++)
        if (NFrom->To == To)
            return NFrom;
    return 0;
}

static int max(const int a, const int b)
{
    return a > b ? a : b;
}
