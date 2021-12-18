#include "libLKH.h"
#include "LKH.h"
#include "Genetic.h"
#include "Heap.h"

/*
 * This file contains the interface functions
 */

static void ResetParameters()
{
    /*
     * ASCENT_CANDIDATES = <integer of at least 2>
     * The number of candidate edges to be associated with each node during the
     * ascent. The candidate set is complemented such that every candidate edge
     * is associated with both its two end nodes.
     */
    AscentCandidates = 50;
    /*
     * BACKBONE_TRIALS = <integer of at least 0>
     * The number of backbone trials in each run.
     */
    BackboneTrials = 0;
    /*
     * BACKTRACKING = { YES | NO }
     * Specifies whether a backtracking k-opt move is to be used as the first
     * move in a sequence of moves (where k = MOVE_TYPE). 
     */
    Backtracking = 0;
    /*
     * CANDIDATE_SET_TYPE = { ALPHA | DELAUNAY [ PURE ] | NEAREST-NEIGHBOR | 
     *                        POPMUSIC | QUADRANT }
     * Specifies the candidate set type.
     * ALPHA is LKH's default type. ALPHA and POPMUSIC are applicable in general.
     * The other types can only be used for instances given by coordinates.
     * The optional suffix PURE for the DELAUNAY type specifies that only
     * edges of the Delaunay graph are used as candidates.
     */
    CandidateSetType = ALPHA;
    DelaunayPure = 0; // in case of DELAUNAY
    Crossover = ERXT;
    DelaunayPartitioning = 0;
    /*
     * EXCESS = <real of at least 0>
     * The maximum alpha-value allowed for any candidate edge is set to
     * EXCESS times the absolute value of the lower bound of a solution
     * tour (determined by the ascent).
     * Default: 1.0/DIMENSION
     */
    Excess = -1;
    /*
     * EXTRA_CANDIDATES = <integer> [ SYMMETRIC ]
     * Number of extra candidate edges to be added to the candidate set
     * of each node. Their candidate set type may be specified after the
     * keyword EXTRA_CANDIDATE_SET_TYPE.
     * The integer may be followed by the keyword SYMMETRIC, signifying
     * that these extra candidate edges is to be complemented such
     * that each of them is associated with both its two end nodes.
     */
    ExtraCandidates = 0;
    ExtraCandidateSetSymmetric = 0;
    /*
     * EXTRA_CANDIDATE_SET_TYPE = { NEAREST-NEIGHBOR | QUADRANT }
     * The candidate set type of extra candidate edges.
     */
    ExtraCandidateSetType = QUADRANT;
    /*
     * GAIN23 = { YES | NO }
     * Specifies whether the Gain23 function is used.
     */
    Gain23Used = 1;
    /*
     * GAIN_CRITERION = { YES | NO }
     * Specifies whether Lin and Kernighan's gain criterion is used.
     */
    GainCriterionUsed = 1;
    /*
     * INITIAL_PERIOD = <integer of at least 0>
     * The length of the first period in the ascent.
     * Default: DIMENSION/2 (but at least 100)
     */
    InitialPeriod = -1;
    /*
     * INITIAL_STEP_SIZE = <integer of at least 1>
     * The initial step size used in the ascent.
     */
    InitialStepSize = 1;
    /*
     * INITIAL_TOUR_ALGORITHM = { BORUVKA | GREEDY | MOORE | NEAREST-NEIGHBOR | 
     *                            QUICK-BORUVKA | SIERPINSKI | WALK }
     * Specifies the algorithm for obtaining an initial tour.
     */
    InitialTourAlgorithm = WALK;
    KarpPartitioning = 0;
    KCenterPartitioning = 0;
    KMeansPartitioning = 0;
    /*
     * KICKS = <integer of at least 0>
     * Specifies the number of times to "kick" a tour found by Lin-Kernighan.
     * Each kick is a random k-swap kick-move. However, if KICKS is zero, then
     * LKH's special kicking strategy, WALK, is used.
     */
    Kicks = 1;
    /*
     * KICK_TYPE = <integer either 0 or at least 4>
     * Specifies the value of k for a random k-swap kick (an extension of the
     * double-bridge move). If KICK_TYPE is zero, then the LKH's special kicking
     * strategy, WALK, is used.
     */
    KickType = 0;
    /*
     * MAX_BREADTH = <integer of at least 0>
     * The maximum number of candidate edges considered at each level of
     * the search for a move.
     */
    MaxBreadth = INT_MAX;
    /*
     * MAX_CANDIDATES = <integer> [ SYMMETRIC ]
     * The maximum number of candidate edges to be associated with each node.
     * The integer may be followed by the keyword SYMMETRIC, signifying
     * that the candidate set is to be complemented such that every candidate
     * edge is associated with both its two end nodes.
     * If MAX_CANDIDATES is zero the candidate sets are made up of the
     * edges represented in the CANDIDATE_FILEs, the INITIAL_TOUR_FILE,
     * the INPUT_TOUR_FILE, the SUBPROBLEM_TOUR_FILE, and the MERGE_TOUR_FILEs.
     */
    MaxCandidates = 5;
    CandidateSetSymmetric = 0;
    /*
     * POPULATION_SIZE = <integer>
     * Specifies the maximum size of the population in the genetic algorithm.
     * Default: 0
     */
    MaxPopulationSize = 0;
    /*
     * MAX_SWAPS = <integer of at least 0>
     * Specifies the maximum number of swaps (flips) allowed in any search
     * for a tour improvement.
     * Default: DIMENSION
     */
    MaxSwaps = -1;
    /*
     * MAX_TRIALS = <integer>
     * The maximum number of trials in each run.
     * Default: DIMENSION
     */
    MaxTrials = -1;
    MoorePartitioning = 0;
    /*
     * MOVE_TYPE = <integer>
     * Specifies the move type to be used as submove in Lin-Kernighan.
     * An integer value k >= 2 signifies that a sequential k-opt move is used.
     */
    MoveType = 5;
    /*
     * NONSEQUENTIAL_MOVE_TYPE = <integer>
     * Specifies the nonsequential move type to be used. A value K >= 4
     * signifies that attempts are made to improve a tour by nonsequential
     * k-opt moves where 4 <= k <= K. Note, however, that the effect depends
     * on the specifications of PATCHING_C and PATCHING_A.
     * Default: (MOVE_TYPE + PATCHING_C + PATCHING_A - 1)
     */
    NonsequentialMoveType = -1;

    /*
     * PATCHING_A = <integer of at least 0> [ RESTRICTED | EXTENDED ]
     * The maximum number of disjoint alternating cycles to be used for
     * patching. An attempt to patch cycles is made if the corresponding
     * non-sequential move is gainful.
     * The integer may be followed by the keyword RESTRICTED or EXTENDED.
     * The keyword RESTRICTED signifies that gainful moves are only
     * considered if all its inclusion edges are candidate edges.
     * The keyword EXTENDED signifies that the non-sequential move need
     * not be gainful if only all its inclusion edges are candidate edges.
     * Default: 1
     */
    PatchingA = 1;
    /*
     * PATCHING_C = <integer of at least 0> [ RESTRICTED | EXTENDED ]
     * The maximum number of disjoint cycles to be patched in an attempt
     * to find a feasible and gainful move. An attempt to patch cycles is
     * made if the corresponding non-sequential move is gainful.
     * The integer may be followed by the keyword RESTRICTED or EXTENDED.
     * The keyword RESTRICTED signifies that gainful moves are only
     * considered if all its inclusion edges are candidate edges.
     * The keyword EXTENDED signifies that the non-sequential move need
     * not be gainful if only all its inclusion edges are candidate edges.
     * Default: 0
     */
    PatchingC = 0;
    PatchingAExtended = 0;
    PatchingARestricted = 0;
    PatchingCExtended = 0;
    PatchingCRestricted = 0;
    /*
     * PRECISION = <integer>
     * The internal precision in the representation of transformed distances:
     *    d[i][j] = PRECISION*c[i][j] + pi[i] + pi[j],
     * where d[i][j], c[i][j], pi[i] and pi[j] are all integral.
     * Default: 100 (which corresponds to 2 decimal places)
     */
    Precision = 100;
    /*
     * POPMUSIC_INITIAL_TOUR = { YES | NO }
     * Specifies whether the best POPMUSIC tour is to be used as intial tour
     * for Lin-Kernighan.
     * Default: NO
     */
    POPMUSIC_InitialTour = 0;
    /*
     * POPMUSIC_MAX_NEIGHBORS = <integer of at least 1>
     * Maximum number of nearest neighbors used as candidates in 3-opt for
     * POPMUSIC.
     * Default: 5
     */
    POPMUSIC_MaxNeighbors = 5;
    /*
     * POPMUSIC_SAMPLE_SIZE = <integer of at least 1>
     * Sample size.
     * Default: 10
     */
    POPMUSIC_SampleSize = 10;
    /*
     * POPMUSIC_SOLUTIONS = <integer of at least 1>
     * Number of solutions to be generated.
     * Default: 50
     */
    POPMUSIC_Solutions = 50;
    /*
     * POPMUSIC_TRIALS = <int>
     * Number of trials used in iterated 3-opt for POPMUSIC.
     * If the value is zero, the number of trials is the size of the subpath
     * to be optimized.
     * Default: 1
     */
    POPMUSIC_Trials = 1;
    /*
     * RECOMBINATION = { IPT | GPX2 }
     * Default: IPT
     */
    Recombination = IPT;
    /*
     * RESTRICTED_SEARCH = { YES | NO }
     * Specifies whether the following search pruning technique is used:
     * The first edge to be broken in a move must not belong to the currently
     * best solution tour. When no solution tour is known, it must not belong
     * to the minimum spanning 1-tree.
     * Default: YES
     */
    RestrictedSearch = 1;
    RohePartitioning = 0;
    /*
     * RUNS = <integer of at least 1>
     * The total number of runs.
     * Default: 10
     */
    Runs = 0;
    /*
     * SEED = <integer>
     * Specifies the initial seed for random number generation. If zero, the
     * seed is derived from the system clock.
     * Default: 1
     */
    Seed = 1;
    SierpinskiPartitioning = 0;
    /*
     * SUBGRADIENT = { YES | NO }
     * Specifies whether the Pi-values should be determined by subgradient
     * optimization.
     * Default: YES
     */
    Subgradient = 1;
    /*
     * SUBSEQUENT_MOVE_TYPE = <integer either 0 or at least 2>
     * Specifies the move type to be used for all moves following the first move
     * in a sequence of moves. The value K >= 2 signifies that a K-opt move is to
     * be used. The value 0 signifies that all moves are of the same type
     * (K = MOVE_TYPE).
     * Default: 0
     */
    SubsequentMoveType = 0;
    /*
     * SUBSEQUENT_PATCHING = { YES | NO }
     * Specifies whether patching is used for moves following the first move
     * in a sequence of moves.
     * Default: YES
     */
    SubsequentPatching = 1;
    /*
     * TIME_LIMIT = <real of at least 0>
     * Specifies a time limit in seconds.
     * Default: DBL_MAX
     */
    TimeLimit = DBL_MAX;
    /*
     * TRACE_LEVEL = <integer>
     * Specifies the level of detail of the output given during the solution
     * process. The value 0 signifies a minimum amount of output. The higher
     * the value is the more information is given.
     * Default: 1
     */
    TraceLevel = 1;

    /* Reset former static local variables */
    sl_Gain23_s1 = 0;
    sl_Gain23_OldReversed = 0;
}

static void AdjustParameters()
{
    int K;

    if (Seed == 0)
        Seed = (unsigned) time(0);
    if (Precision == 0)
        Precision = 100;
    if (InitialStepSize == 0)
        InitialStepSize = 1;
    if (MaxSwaps < 0)
        MaxSwaps = Dimension;
    if (KickType > Dimension / 2)
        KickType = Dimension / 2;
    if (Runs == 0)
        Runs = 10;
    if (MaxCandidates > Dimension - 1)
        MaxCandidates = Dimension - 1;
    if (ExtraCandidates > Dimension - 1)
        ExtraCandidates = Dimension - 1;
    if (AscentCandidates > Dimension - 1)
        AscentCandidates = Dimension - 1;
    if (InitialPeriod < 0) {
        InitialPeriod = Dimension / 2;
        if (InitialPeriod < 100)
            InitialPeriod = 100;
    }
    if (Excess < 0)
        Excess = 1.0 / Dimension;
    if (MaxTrials == -1)
        MaxTrials = Dimension;
    MakeHeap(Dimension);
    if (POPMUSIC_MaxNeighbors > Dimension - 1)
        POPMUSIC_MaxNeighbors = Dimension - 1;
    if (POPMUSIC_SampleSize > Dimension)
        POPMUSIC_SampleSize = Dimension;
    if (SubsequentMoveType == 0)
        SubsequentMoveType = MoveType;
    K = MoveType >= SubsequentMoveType
        || !SubsequentPatching ? MoveType : SubsequentMoveType;
    if (PatchingC > K)
        PatchingC = K;
    if (PatchingA > 1 && PatchingA >= PatchingC)
        PatchingA = PatchingC > 2 ? PatchingC - 1 : 1;
    if (NonsequentialMoveType == -1 ||
            NonsequentialMoveType > K + PatchingC + PatchingA - 1)
        NonsequentialMoveType = K + PatchingC + PatchingA - 1;
    if (PatchingC >= 1) {
        BestMove = BestSubsequentMove = BestKOptMove;
        if (!SubsequentPatching && SubsequentMoveType <= 5) {
            MoveFunction BestOptMove[] =
            { 0, 0, Best2OptMove, Best3OptMove,
                Best4OptMove, Best5OptMove
            };
            BestSubsequentMove = BestOptMove[SubsequentMoveType];
        }
    } else {
        MoveFunction BestOptMove[] = { 0, 0, Best2OptMove, Best3OptMove,
            Best4OptMove, Best5OptMove
        };
        BestMove = MoveType <= 5 ? BestOptMove[MoveType] : BestKOptMove;
        BestSubsequentMove = SubsequentMoveType <= 5 ?
            BestOptMove[SubsequentMoveType] : BestKOptMove;
    }
}

static void ReadCoords(struct NodeCoords const * coords)
{
    Node *Prev = 0, *N = 0;
    int i;

    NodeSet = (Node *) calloc(Dimension + 1, sizeof(Node));
    for (i = 1; i <= Dimension; i++, Prev = N) {
        N = &NodeSet[i];
        N->V = 1;
        N->X = coords[i-1].x;
        N->Y = coords[i-1].y;

        if (i == 1)
            FirstNode = N;
        else
            Link(Prev, N);
        N->Id = i;
    }
    Link(N, FirstNode);
}

int const *run(int dimension, struct NodeCoords const * coords)
{
    GainType Cost;
    double Time, LastTime;

    ResetParameters();

    Dimension = dimension;
    if (Dimension < 3)
        eprintf("DIMENSION < 3 or not specified");

    /* Read the specification of the problem */
    StartTime = LastTime = GetTime();
    MaxMatrixDimension = 20000;
    MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
        MergeWithTourGPX2;
    FreeStructures();
    FirstNode = 0;
    ReadCoords(coords);
    Swaps = 0;
    AdjustParameters();

    AllocateStructures();
    CreateCandidateSet();
    InitializeStatistics();

    if (Norm != 0)
        BestCost = PLUS_INFINITY;
    else {
        /* The ascent has solved the problem! */
        BestCost = (GainType) LowerBound;
        UpdateStatistics(BestCost, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        Runs = 0;
    }

    /* Find a specified number (Runs) of local optima */
    for (Run = 1; Run <= Runs; Run++) {
        LastTime = GetTime();
        if (LastTime - StartTime >= TimeLimit) {
            if (TraceLevel >= 1)
                printff("*** Time limit exceeded ***\n");
            break;
        }
        Cost = FindTour();      /* using the Lin-Kernighan heuristic */
        if (MaxPopulationSize > 1) {
            /* Genetic algorithm */
            int i;
            for (i = 0; i < PopulationSize; i++) {
                GainType OldCost = Cost;
                Cost = MergeTourWithIndividual(i);
                if (TraceLevel >= 1 && Cost < OldCost) {
                    printff("  Merged with %d: Cost = " GainFormat "\n", i + 1,
                            Cost);
                }
            }
            if (!HasFitness(Cost)) {
                if (PopulationSize < MaxPopulationSize) {
                    AddToPopulation(Cost);
                    if (TraceLevel >= 1)
                        PrintPopulation();
                } else if (Cost < Fitness[PopulationSize - 1]) {
                    i = ReplacementIndividual(Cost);
                    ReplaceIndividualWithTour(i, Cost);
                    if (TraceLevel >= 1)
                        PrintPopulation();
                }
            }
        } else if (Run > 1)
            Cost = MergeTourWithBestTour();
        if (Cost < BestCost) {
            BestCost = Cost;
            RecordBetterTour();
            RecordBestTour();
        }
        Time = fabs(GetTime() - LastTime);
        UpdateStatistics(Cost, Time);
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
            printff("Run %d: Cost = " GainFormat, Run, Cost);
            printff(", Time = %0.2f sec.\n\n", Time);
        }
        if (PopulationSize >= 2 &&
            (PopulationSize == MaxPopulationSize ||
             Run >= 2 * MaxPopulationSize) && Run < Runs) {
            Node *N;
            int Parent1, Parent2;
            Parent1 = LinearSelection(PopulationSize, 1.25);
            do
                Parent2 = LinearSelection(PopulationSize, 1.25);
            while (Parent2 == Parent1);
            ApplyCrossover(Parent1, Parent2);
            N = FirstNode;
            do {
                int d = C(N, N->Suc);
                AddCandidate(N, N->Suc, d, INT_MAX);
                AddCandidate(N->Suc, N, d, INT_MAX);
                N = N->InitialSuc = N->Suc;
            }
            while (N != FirstNode);
        }
        SRandom(++Seed);
    }
    PrintStatistics();
    return BestTour;
}
