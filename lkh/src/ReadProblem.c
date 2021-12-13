#include "LKH.h"
#include "Heap.h"

/*
 * The ReadProblem function reads the problem data in TSPLIB format from the
 * file specified in the parameter file (PROBLEM_FILE).
 *
 * The following description of the file format is extracted from the TSPLIB
 * documentation.
 *
 * The file consists of a specification part and a data part. The specification
 * part contains information on the file format and on its contents. The data
 * part contains explicit data.
 *
 * (1) The specification part
 *
 * All entries in this section are of the form <keyword> : <value>, where
 * <keyword> denotes an alphanumerical keyword and <value> denotes
 * alphanumerical or numerical data. The terms <string>, <integer> and <real>
 * denote character string, integer or real data, respectively. The order of
 * specification of the keywords in the data file is arbitrary (in principle),
 * but must be consistent, i.e., whenever a keyword is specified, all
 * necessary information for the correct interpretation of the keyword has to
 * be known.
 *
 * Below is given a list of all available keywords.
 *
 * EDGE_DATA_FORMAT : <string>
 * Describes the format in which the edges of a graph are given, if the
 * graph is not complete. The values are
 * EDGE_LIST    The graph is given by an edge list
 * ADJ_LIST     The graph is given by an adjacency list
 *
 * GRID_SIZE : <real>
 * The grid size for toroidal instances.
 * Default: 1000000.0
 *
 * EOF
 * Terminates input data. The entry is optional.
 *
 * (2) The data part
 *
 * Depending on the choice of specifications some additional data may be
 * required. These data are given corresponding data sections following the
 * specification part. Each data section begins with the corresponding
 * keyword. The length of the sectionis either explicitly known form the
 * format specification, or the section is terminated by an appropriate
 * end-of-section identifier.
 *
 * NODE_COORD_SECTION :
 * Node coordinates are given in this section. Each line is of the form
 *
 *      <integer> <real> <real>
 *
 * The real numbers are the associated coordinates.
 */

static const char Delimiters[] = " :=\n\t\r\f\v\xef\xbb\xbf";
static void CreateNodes(void);
static void Read_GRID_SIZE(void);
static void Read_NODE_COORD_SECTION(void);

static FILE *ProblemFile;

void ReadProblem()
{
    int i, K;
    char *Line, *Keyword;

    ProblemFile = 0;

    if (!(ProblemFile = fopen("input", "r")))
        eprintf("Cannot open PROBLEM_FILE: \"input\"");

    FreeStructures();
    FirstNode = 0;
    GridSize = 1000000.0;
    while ((Line = ReadLine(ProblemFile))) {
        if (!(Keyword = strtok(Line, Delimiters)))
            continue;
        for (i = 0; i < (int) strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
        if (!strcmp(Keyword, "EOF"))
            break;
        else if (!strcmp(Keyword, "GRID_SIZE"))
            Read_GRID_SIZE();
        else if (!strcmp(Keyword, "NODE_COORD_SECTION"))
            Read_NODE_COORD_SECTION();
        else
            eprintf("Unknown keyword: %s", Keyword);
    }
    Swaps = 0;

    /* Adjust parameters */
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
    if (TraceLevel >= 1) {
        printff("done\n");
    }
    fclose(ProblemFile);
    free(LastLine);
    LastLine = 0;
}

static void CreateNodes()
{
    Node *Prev = 0, *N = 0;
    int i;

    NodeSet = (Node *) calloc(Dimension + 1, sizeof(Node));
    for (i = 1; i <= Dimension; i++, Prev = N) {
        N = &NodeSet[i];
        if (i == 1)
            FirstNode = N;
        else
            Link(Prev, N);
        N->Id = i;
    }
    Link(N, FirstNode);
}

static void Read_GRID_SIZE()
{
    char *Token = strtok(0, Delimiters);

    if (!Token || !sscanf(Token, "%lf", &GridSize))
        eprintf("GRID_SIZE: real expected");
    if (GridSize < 0)
        eprintf("GRID_SIZE: non-negative real expected");
}

static void Read_NODE_COORD_SECTION()
{
    Node *N;
    int Id, i;

    if (!FirstNode)
        CreateNodes();
    N = FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != FirstNode);
    for (i = 1; i <= Dimension; i++) {
        if (!fscanint(ProblemFile, &Id))
            eprintf("NODE_COORD_SECTION: Missing nodes");
        if (Id <= 0 || Id > Dimension)
            eprintf("NODE_COORD_SECTION: Node number out of range: %d",
                    Id);
        N = &NodeSet[Id];
        if (N->V == 1)
            eprintf("NODE_COORD_SECTION: Node number occurs twice: %d",
                    N->Id);
        N->V = 1;
        if (!fscanf(ProblemFile, "%lf", &N->X))
            eprintf("NODE_COORD_SECTION : Missing X-coordinate");
        if (!fscanf(ProblemFile, "%lf", &N->Y))
            eprintf("NODE_COORD_SECTION: Missing Y-coordinate");
    }
    N = FirstNode;
    do
        if (!N->V && N->Id <= Dimension)
            break;
    while ((N = N->Suc) != FirstNode);
    if (!N->V)
        eprintf("NODE_COORD_SECTION: No coordinates given for node %d",
                N->Id);
}
