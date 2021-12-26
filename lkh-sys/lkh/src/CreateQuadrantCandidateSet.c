#include "LKH.h"

/*
 * This file contains the CreateQuadrantCandidateSet function
 * and the CreateNearestNeighborCandidateSet function.
 */

static void NearestQuadrantNeighbors(Node * N, int Q, int K);
static int Contains(Node * T, int Q, Node * N);
static int BoxOverlaps(Node * T, int Q, Node * N);
static void ComputeBounds(int start, int end);

typedef int (*ContainsFunction) (Node * T, int Q, Node * N);
typedef int (*BoxOverlapsFunction) (Node * T, int Q, Node * N);

static Node **KDTree;
static Candidate *CandidateSet;
static double *XMin, *XMax, *YMin, *YMax;
static int Candidates, Radius;
static int Level = 0;

/*
 * The CreateQuadrantCandidateSet function creates for each node 
 * a candidate set consisting of the K/L least costly neighbor edges
 * in each of the L geometric quadrants around the node, where L 
 * is 4 for 2-D instances, and 8 for 3-D instances. 
 * If these totals less than K nodes, the candidate set is 
 * augmented by the nearest remaining nodes overall to bring 
 * the total up to K.
 *
 * The function is called from the CreateCandidateSet and
 * AddExtraCandidates functions.
 */

void CreateQuadrantCandidateSet(int K)
{
    Node *From, *To;
    int L, Q, CandPerQ, Added, i;

    if (K <= 0)
        return;
    if (TraceLevel >= 2)
        printff("Creating quadrant candidate set ... ");
    KDTree = BuildKDTree(1);
    XMin = (double *) malloc((1 + Dimension) * sizeof(double));
    XMax = (double *) malloc((1 + Dimension) * sizeof(double));
    YMin = (double *) malloc((1 + Dimension) * sizeof(double));
    YMax = (double *) malloc((1 + Dimension) * sizeof(double));
    ComputeBounds(0, Dimension - 1);
    L = 4;
    CandPerQ = K / L;
    CandidateSet = (Candidate *) malloc((K + 1) * sizeof(Candidate));

    From = FirstNode;
    do {
        if (FixedOrCommonCandidates(From) == 2)
            continue;
        Added = 0;
        for (Q = 1; Q <= L; Q++) {
            NearestQuadrantNeighbors(From, Q, CandPerQ);
            for (i = 0; i < Candidates; i++) {
                To = CandidateSet[i].To;
                if (AddCandidate(From, To, D(From, To), 1))
                    Added++;
            }
        }
        if (K > Added) {
            NearestQuadrantNeighbors(From, 0, K - Added);
            for (i = 0; i < Candidates; i++) {
                To = CandidateSet[i].To;
                AddCandidate(From, To, D(From, To), 2);
            }
        }
    } while ((From = From->Suc) != FirstNode);

    free(CandidateSet);
    free(KDTree);
    free(XMin);
    free(XMax);
    free(YMin);
    free(YMax);
    if (Level == 0) {
        ResetCandidateSet();
        AddTourCandidates();
        if (CandidateSetSymmetric)
            SymmetrizeCandidateSet();
        if (TraceLevel >= 2)
            printff("done\n");
    }
}

/*
 * The CreateNearestNeighborCandidateSet function creates for each node
 * a candidate set consisting of the K least costly neighbor edges.
 *
 * The function is called from the CreateCandidateSet and 
 * AddExtraCandidates functions.
 */

void CreateNearestNeighborCandidateSet(int K)
{
    Node *From, *To;
    int i;

    if (TraceLevel >= 2)
        printff("Creating nearest neighbor candidate set ... ");
    KDTree = BuildKDTree(1);
    XMin = (double *) malloc((1 + Dimension) * sizeof(double));
    XMax = (double *) malloc((1 + Dimension) * sizeof(double));
    YMin = (double *) malloc((1 + Dimension) * sizeof(double));
    YMax = (double *) malloc((1 + Dimension) * sizeof(double));
    ComputeBounds(0, Dimension - 1);
    CandidateSet = (Candidate *) malloc((K + 1) * sizeof(Candidate));

    From = FirstNode;
    do {
        NearestQuadrantNeighbors(From, 0, K);
        for (i = 0; i < Candidates; i++) {
            To = CandidateSet[i].To;
            AddCandidate(From, To, D(From, To), 1);
        }
    } while ((From = From->Suc) != FirstNode);

    free(CandidateSet);
    free(KDTree);
    free(XMin);
    free(XMax);
    free(YMin);
    free(YMax);
    if (Level == 0) {
        ResetCandidateSet();
        AddTourCandidates();
        if (CandidateSetSymmetric)
            SymmetrizeCandidateSet();
        if (TraceLevel >= 2)
            printff("done\n");
    }
}

/*
 * The ComputeBounds function computes the bounding boxes
 * for the K-d tree nodes.
 */

static void ComputeBounds(int start, int end)
{
    if (start <= end) {
        int mid = (start + end) / 2, i;
        Node *T = KDTree[mid];
        XMin[T->Id] = YMin[T->Id] = DBL_MAX;
        XMax[T->Id] = YMax[T->Id] = -DBL_MAX;
        for (i = start; i <= end; i++) {
            Node *N = KDTree[i];
            if (N == T)
                continue;
            if (N->X < XMin[T->Id])
                XMin[T->Id] = N->X;
            if (N->X > XMax[T->Id])
                XMax[T->Id] = N->X;
            if (N->Y < YMin[T->Id])
                YMin[T->Id] = N->Y;
            if (N->Y > YMax[T->Id])
                YMax[T->Id] = N->Y;
        }
        ComputeBounds(start, mid - 1);
        ComputeBounds(mid + 1, end);
    }
}


#define Coord(N, axis) (axis == 0 ? (N)->X : (N)->Y)

/*
 * The Contains2D function returns 1 if T belongs to 2-D quadrant Q
 * relative to N; otherwise 0.
 *
 *          Q = 2 | Q = 1
 *          ===== N =====
 *          Q = 3 | Q = 4
 *
 */

static int Contains(Node * T, int Q, Node * N)
{
    switch (Q) {
    case 1:
        return T->X >= N->X && T->Y >= N->Y;
    case 2:
        return T->X <= N->X && T->Y >= N->Y;
    case 3:
        return T->X <= N->X && T->Y <= N->Y;
    case 4:
        return T->X >= N->X && T->Y <= N->Y;
    default:
        return 1;
    }
}

/*
 * The BoxOverlaps2D function returns 1 if T's bounding box 
 * overlaps the 2-D quadrant Q relative to N; otherwise 0.
 */

static int BoxOverlaps(Node * T, int Q, Node * N)
{
    switch (Q) {
    case 1:
        return XMax[T->Id] >= N->X && YMax[T->Id] >= N->Y;
    case 2:
        return XMin[T->Id] <= N->X && YMax[T->Id] >= N->Y;
    case 3:
        return XMin[T->Id] <= N->X && YMin[T->Id] <= N->Y;
    case 4:
        return XMax[T->Id] >= N->X && YMin[T->Id] <= N->Y;
    default:
        return 1;
    }
}


/*
 * The Overlaps function returns 1 if High is zero and the half 
 * plane to the left of a point T overlaps quadrant Q relative 
 * to a point N.
 * If High is not zero the Overlaps function returns 1 if the 
 * half plane to the right of T overlaps quadrant Q relative to N.
 * Otherwise the function returns 0.
 *
 * The directions are relative to the given coordinate axis.
 * The parameter diff is <= 0 if T is to the left of N, 
 * and >= 0 if T is to the right of N.   
 */

static int Overlaps(int Q, double diff, int High, int axis)
{
    switch (Q) {
    case 1:
        return High || diff >= 0;
    case 2:
        return axis == 0 ? !High || diff <= 0 : High || diff >= 0;
    case 3:
        return axis <= 1 ? !High || diff <= 0 : High || diff >= 0;
    case 4:
        return axis == 1 ? !High || diff <= 0 : High || diff >= 0;
    case 5:
        return axis <= 1 ? High || diff >= 0 : !High || diff <= 0;
    case 6:
        return axis == 1 ? High || diff >= 0 : !High || diff <= 0;
    case 7:
        return !High || diff <= 0;
    case 8:
        return axis == 0 ? High || diff >= 0 : !High || diff <= 0;
    default:
        return 1;
    }
}

static int InCandidateSet(Node * N, Node * T)
{
    int i;
    for (i = 0; i < Candidates; i++)
        if (CandidateSet[i].To == T)
            return 1;
    return IsCandidate(N, T);
}

/*
 * The NQN function searches KDTree[start:end] in an attempt to
 * find the K quad-nearest neighbors in quadrant Q relative to N. 
 *
 * The function is called from the NearestQuadrantNeighbors 
 * and the CreateNearestNeighborCandidateSet function.
 */

static void NQN(Node * N, int Q, int start, int end, int K)
{
    int mid = (start + end) / 2, d;
    Node *T = KDTree[mid], P = { 0 };
    int axis = T->Axis;

    if (K == 0 || N->FixedTo2)
        return;
    if (start <= end && T != N && !T->FixedTo2 &&
        IsPossibleCandidate(N, T) &&
        Contains(T, Q, N) &&
        !InCandidateSet(N, T) &&
        (c(N, T) - N->Pi - T->Pi <= Radius) &&
        (d = Distance(N, T) * Precision) <= Radius) {
        int i = Candidates;
        while (--i >= 0 && d < CandidateSet[i].Cost)
            CandidateSet[i + 1] = CandidateSet[i];
        CandidateSet[i + 1].To = T;
        CandidateSet[i + 1].Cost = d;
        if (Candidates < K)
            Candidates++;
        if (Candidates == K)
            Radius = CandidateSet[Candidates - 1].Cost;
    }
    if (start < end && BoxOverlaps(T, Q, N)) {
        double diff = Coord(T, axis) - Coord(N, axis);
        P.X = axis == 0 ? T->X : N->X;
        P.Y = axis == 1 ? T->Y : N->Y;
        P.Pi = 0;
        if (diff >= 0) {
            if (Overlaps(Q, diff, 0, axis))
                NQN(N, Q, start, mid - 1, K);
            if (Overlaps(Q, diff, 1, axis) &&
                (c(N, &P) - N->Pi <= Radius) &&
                Distance(N, &P) * Precision <= Radius)
                NQN(N, Q, mid + 1, end, K);
        } else {
            if (Overlaps(Q, diff, 1, axis))
                NQN(N, Q, mid + 1, end, K);
            if (Overlaps(Q, diff, 0, axis) &&
                (c(N, &P) - N->Pi <= Radius) &&
                Distance(N, &P) * Precision <= Radius)
                NQN(N, Q, start, mid - 1, K);
        }
    }
}

/*
 * The NearestQuadrantNeighbors function searches the K-d tree 
 * in an attempt to find the K quad-nearest neighbors in 
 * quadrant Q relative to N.
 * If Q = 0, the funtion computes the K nearest neighbors to N.
 */

static void NearestQuadrantNeighbors(Node * N, int Q, int K)
{
    Candidates = 0;
    Radius = INT_MAX;
    NQN(N, Q, 0, Dimension - 1, K);
}
