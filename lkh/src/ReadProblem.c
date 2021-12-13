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
 * EDGE-WEIGHT_FORMAT : <string>
 * Describes the format of the edge weights if they are given explicitly.
 * The values are
 * FUNCTION         Weights are given by a function (see above)
 * FULL_MATRIX      Weights are given by a full matrix
 * UPPER_ROW        Upper triangular matrix
 *                      (row-wise without diagonal entries)
 * LOWER_ROW        Lower triangular matrix
 *                      (row-wise without diagonal entries)
 * UPPER_DIAG_ROW   Upper triangular matrix
 *                      (row-wise including diagonal entries)
 * LOWER_DIAG_ROW   Lower triangular matrix
 *                      (row-wise including diagonal entries)
 * UPPER_COL        Upper triangular matrix
 *                      (column-wise without diagonal entries)
 * LOWER_COL        Lower triangular matrix
 *                      (column-wise without diagonal entries)
 * UPPER_DIAG_COL   Upper triangular matrix
 *                      (column-wise including diagonal entries)
 * LOWER_DIAG_COL   Lower triangular matrix
 *                      (column-wise including diagonal entries)
 *
 * EDGE_DATA_FORMAT : <string>
 * Describes the format in which the edges of a graph are given, if the
 * graph is not complete. The values are
 * EDGE_LIST    The graph is given by an edge list
 * ADJ_LIST     The graph is given by an adjacency list
 *
 * DISPLAY_DATA_TYPE : <string>
 * Specifies how a graphical display of the nodes can be obtained.
 * The values are
 * COORD_DISPLAY    Display is generated from the node coordinates
 * TWOD_DISPLAY     Explicit coordinates in 2-D are given
 * NO_DISPLAY       No graphical display is possible
 *
 * The default value is COORD_DISPLAY if node coordinates are specifies and
 * NO_DISPLAY otherwise. In the current implementation, however, the value
 * has no significance.
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
 *
 * EDGE_DATA_SECTION :
 * Edges of the graph are specified in either of the two formats allowed in
 * the EDGE_DATA_FORMAT entry. If the type is EDGE_LIST, then the edges are
 * given as a sequence of lines of one of the forms
 *
 *      <integer> <integer>
 *      <integer> <integer> <integer>
 *
 * each entry giving the terminal nodes of some edge, and if three integers are
 * given, the last one specifies its weight. The list is terminated by a -1.
 * If the type is ADJ_LIST, the section consists of adjacency lists for nodes.
 * The adjacency list of a node x is specified as
 *
 *      <integer> <integer> ... <integer> -1
 *
 * where the first integer gives the number of node x and the following
 * integers (terminated by -1) the numbers of the nodes adjacent to x.
 * The list of adjacency lists are terminated by an additional -1.
 *
 * FIXED_EDGES_SECTION :
 * In this section, edges are listed that are required to appear in each
 * solution to the problem. The edges to be fixed are given in the form
 * (per line)
 *
 *      <integer> <integer>
 *
 * meaning that the edge (arc) from the first node to the second node has
 * to be contained in a solution. This section is terminated by a -1.
 *
 * DISPLAY_DATA_SECTION :
 * If DISPLAY_DATA_TYPE is TWOD_DISPLAY, the 2-dimensional coordinates from
 * which a display can be generated are given in the form (per line)
 *
 *      <integer> <real> <real>
 *
 * The integers specify the respective nodes and the real numbers give the
 * associated coordinates. The contents of this section, however, has no
 * significance in the current implementation.
 *
 * EDGE_WEIGHT_SECTION :
 * The edge weights are given in the format specifies by the EDGE_WEIGHT_FORMAT
 * entry. At present, all explicit data are integral and is given in one of the
 * (self-explanatory) matrix formats, with explicitly known lengths.
 */

static const char Delimiters[] = " :=\n\t\r\f\v\xef\xbb\xbf";
static void CheckSpecificationPart(void);
static char *Copy(char *S);
static void CreateNodes(void);
static int FixEdge(Node * Na, Node * Nb);
static void Read_DISPLAY_DATA_SECTION(void);
static void Read_DISPLAY_DATA_TYPE(void);
static void Read_EDGE_DATA_FORMAT(void);
static void Read_EDGE_WEIGHT_FORMAT(void);
static void Read_EDGE_WEIGHT_SECTION(void);
static void Read_FIXED_EDGES_SECTION(void);
static void Read_GRID_SIZE(void);
static void Read_NODE_COORD_SECTION(void);
static int TwoDWeightType(void);
static int ThreeDWeightType(void);

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
    WeightFormat = -1;
    Type = EdgeWeightType = EdgeWeightFormat = 0;
    EdgeDataFormat = NodeCoordType = DisplayDataType = 0;
    GridSize = 1000000.0;
    while ((Line = ReadLine(ProblemFile))) {
        if (!(Keyword = strtok(Line, Delimiters)))
            continue;
        for (i = 0; i < (int) strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
        if (!strcmp(Keyword, "DISPLAY_DATA_SECTION"))
            Read_DISPLAY_DATA_SECTION();
        else if (!strcmp(Keyword, "DISPLAY_DATA_TYPE"))
            Read_DISPLAY_DATA_TYPE();
        else if (!strcmp(Keyword, "EDGE_DATA_FORMAT"))
            Read_EDGE_DATA_FORMAT();
        else if (!strcmp(Keyword, "EDGE_WEIGHT_FORMAT"))
            Read_EDGE_WEIGHT_FORMAT();
        else if (!strcmp(Keyword, "EDGE_WEIGHT_SECTION"))
            Read_EDGE_WEIGHT_SECTION();
        else if (!strcmp(Keyword, "EOF"))
            break;
        else if (!strcmp(Keyword, "FIXED_EDGES_SECTION"))
            Read_FIXED_EDGES_SECTION();
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

static int TwoDWeightType()
{
    return 1;
}

static int ThreeDWeightType()
{
    return 0;
}

static void CheckSpecificationPart()
{
    if (WeightFormat != -1
        && WeightFormat != FUNCTION)
        eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if (CandidateSetType == DELAUNAY && !TwoDWeightType()
        && MaxCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = DELAUNAY");
    if (CandidateSetType == NN && !TwoDWeightType()
        && !ThreeDWeightType() && MaxCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = "
             "NEAREST-NEIGHBOR");
    if (CandidateSetType == QUADRANT && !TwoDWeightType()
        && !ThreeDWeightType() && MaxCandidates + ExtraCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = QUADRANT");
    if (ExtraCandidateSetType == NN && !TwoDWeightType()
        && !ThreeDWeightType() && ExtraCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for EXTRA_CANDIDATE_SET_TYPE = "
             "NEAREST-NEIGHBOR");
    if (ExtraCandidateSetType == QUADRANT && !TwoDWeightType()
        && !ThreeDWeightType()
        && ExtraCandidates > 0)
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for EXTRA_CANDIDATE_SET_TYPE = "
             "QUADRANT");
    if (InitialTourAlgorithm == QUICK_BORUVKA && !TwoDWeightType()
        && !ThreeDWeightType())
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "QUICK-BORUVKA");
    if (InitialTourAlgorithm == SIERPINSKI && !TwoDWeightType())
        eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "SIERPINSKI");
    if (DelaunayPartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for DELAUNAY specification");
    if (KarpPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for KARP specification");
    if (KCenterPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for K-CENTER specification");
    if (KMeansPartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for K-MEANS specification");
    if (MoorePartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for MOORE specification");
    if (RohePartitioning && !TwoDWeightType() && !ThreeDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for ROHE specification");
    if (SierpinskiPartitioning && !TwoDWeightType())
        eprintf("Illegal EDGE_WEIGHT_TYPE for SIERPINSKI specification");
}

static char *Copy(char *S)
{
    char *Buffer;

    if (!S || strlen(S) == 0)
        return 0;
    Buffer = (char *) malloc(strlen(S) + 1);
    strcpy(Buffer, S);
    return Buffer;
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

static int FixEdge(Node * Na, Node * Nb)
{
    if (!Na->FixedTo1 || Na->FixedTo1 == Nb)
        Na->FixedTo1 = Nb;
    else if (!Na->FixedTo2 || Na->FixedTo2 == Nb)
        Na->FixedTo2 = Nb;
    else
        return 0;
    if (!Nb->FixedTo1 || Nb->FixedTo1 == Na)
        Nb->FixedTo1 = Na;
    else if (!Nb->FixedTo2 || Nb->FixedTo1 == Na)
        Nb->FixedTo2 = Na;
    else
        return 0;
    return 1;
}

static void Read_DISPLAY_DATA_SECTION()
{
    Node *N;
    int Id, i;

    CheckSpecificationPart();
    if (!DisplayDataType || strcmp(DisplayDataType, "TWOD_DISPLAY"))
        eprintf
            ("DISPLAY_DATA_SECTION conflicts with DISPLAY_DATA_TYPE: %s",
             DisplayDataType);
    if (!FirstNode)
        CreateNodes();
    N = FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != FirstNode);
    for (i = 1; i <= Dimension; i++) {
        if (!fscanint(ProblemFile, &Id))
            eprintf("DIPLAY_DATA_SECTION: Missing nodes");
        if (Id <= 0 || Id > Dimension)
            eprintf("DIPLAY_DATA_SECTION: Node number out of range: %d",
                    Id);
        N = &NodeSet[Id];
        if (N->V == 1)
            eprintf("DIPLAY_DATA_SECTION: Node number occurs twice: %d",
                    N->Id);
        N->V = 1;
        if (!fscanf(ProblemFile, "%lf", &N->X))
            eprintf("DISPLAY_DATA_SECTION: Missing X-coordinate");
        if (!fscanf(ProblemFile, "%lf", &N->Y))
            eprintf("DISPLAY_DATA_SECTION: Missing Y-coordinate");
    }
    N = FirstNode;
    do
        if (!N->V)
            break;
    while ((N = N->Suc) != FirstNode);
    if (!N->V)
        eprintf("DIPLAY_DATA_SECTION: No coordinates given for node %d",
                N->Id);
}

static void Read_DISPLAY_DATA_TYPE()
{
    unsigned int i;

    if (!(DisplayDataType = Copy(strtok(0, Delimiters))))
        eprintf("DISPLAY_DATA_TYPE: string expected");
    for (i = 0; i < strlen(DisplayDataType); i++)
        DisplayDataType[i] = (char) toupper(DisplayDataType[i]);
    if (strcmp(DisplayDataType, "COORD_DISPLAY") &&
        strcmp(DisplayDataType, "TWOD_DISPLAY") &&
        strcmp(DisplayDataType, "NO_DISPLAY"))
        eprintf("Unknown DISPLAY_DATA_TYPE: %s", DisplayDataType);
}

static void Read_EDGE_DATA_FORMAT()
{
    unsigned int i;

    if (!(EdgeDataFormat = Copy(strtok(0, Delimiters))))
        eprintf("EDGE_DATA_FORMAT: string expected");
    for (i = 0; i < strlen(EdgeDataFormat); i++)
        EdgeDataFormat[i] = (char) toupper(EdgeDataFormat[i]);
    if (strcmp(EdgeDataFormat, "EDGE_LIST") &&
        strcmp(EdgeDataFormat, "ADJ_LIST"))
        eprintf("Unknown EDGE_DATA_FORMAT: %s", EdgeDataFormat);
}

/*
static void Read_EDGE_DATA_SECTION()
{
    Node *Ni, *Nj;
    int i, j, W = 0, WithWeights = 0, FirstLine = 1;
    char *Line;

    CheckSpecificationPart();
    if (!EdgeDataFormat)
        eprintf("Missing EDGE_DATA_FORMAT specification");
    if (!FirstNode)
        CreateNodes();
    if (!strcmp(EdgeDataFormat, "EDGE_LIST")) {
        Line = ReadLine(ProblemFile);
        if (sscanf(Line, "%d %d %d\n", &i, &j, &W) == 3)
            WithWeights = 1;
        while (i != -1) {
            if (i <= 0 ||
                i > Dimension)
                eprintf("EDGE_DATA_SECTION: Node number out of range: %d",
                        i);
            if (!FirstLine)
                fscanint(ProblemFile, &j);
            if (j <= 0 || j > Dimension)
                eprintf
                    ("EDGE_DATA_SECTION: Node number out of range: %d",
                     j);
            if (i == j)
                eprintf("EDGE_DATA_SECTION: Illegal edge: %d to %d",
                        i, j);
            Ni = &NodeSet[i];
            Nj = &NodeSet[j];
            if (WithWeights) {
                if (!FirstLine)
                    fscanint(ProblemFile, &W);
                W *= Precision;
            }
            AddCandidate(Ni, Nj, W, 1);
            AddCandidate(Nj, Ni, W, 1);
            FirstLine = 0;
            if (!fscanint(ProblemFile, &i))
                i = -1;
        }
    } else if (!strcmp(EdgeDataFormat, "ADJ_LIST")) {
        if (!fscanint(ProblemFile, &i))
            i = -1;
        while (i != -1) {
            if (i <= 0 || Dimension)
                eprintf
                    ("EDGE_DATA_SECTION: Node number out of range: %d",
                     i);
            Ni = &NodeSet[i];
            fscanint(ProblemFile, &j);
            while (j != -1) {
                if (j <= 0 || Dimension)
                    eprintf
                        ("EDGE_DATA_SECTION: Node number out of range: %d",
                         j);
                if (i == j)
                    eprintf("EDGE_DATA_SECTION: Illgal edge: %d to %d",
                            i, j);
                Nj = &NodeSet[j];
                AddCandidate(Ni, Nj, 0, 1);
                AddCandidate(Nj, Ni, 0, 1);
                fscanint(ProblemFile, &j);
            }
            fscanint(ProblemFile, &i);
        }
    } else
        eprintf("EDGE_DATA_SECTION: No EDGE_DATA_FORMAT specified");
    WeightType = -1;
    MaxCandidates = ExtraCandidates = 0;
    Distance = WithWeights ? Distance_LARGE : Distance_1;
}
*/

static void Read_EDGE_WEIGHT_FORMAT()
{
    unsigned int i;

    if (!(EdgeWeightFormat = Copy(strtok(0, Delimiters))))
        eprintf("EDGE_WEIGHT_FORMAT: string expected");
    for (i = 0; i < strlen(EdgeWeightFormat); i++)
        EdgeWeightFormat[i] = (char) toupper(EdgeWeightFormat[i]);
    if (!strcmp(EdgeWeightFormat, "FUNCTION"))
        WeightFormat = FUNCTION;
    else if (!strcmp(EdgeWeightFormat, "FULL_MATRIX"))
        WeightFormat = FULL_MATRIX;
    else if (!strcmp(EdgeWeightFormat, "UPPER_ROW"))
        WeightFormat = UPPER_ROW;
    else if (!strcmp(EdgeWeightFormat, "LOWER_ROW"))
        WeightFormat = LOWER_ROW;
    else if (!strcmp(EdgeWeightFormat, "UPPER_DIAG_ROW"))
        WeightFormat = UPPER_DIAG_ROW;
    else if (!strcmp(EdgeWeightFormat, "LOWER_DIAG_ROW"))
        WeightFormat = LOWER_DIAG_ROW;
    else if (!strcmp(EdgeWeightFormat, "UPPER_COL"))
        WeightFormat = UPPER_COL;
    else if (!strcmp(EdgeWeightFormat, "LOWER_COL"))
        WeightFormat = LOWER_COL;
    else if (!strcmp(EdgeWeightFormat, "UPPER_DIAG_COL"))
        WeightFormat = UPPER_DIAG_COL;
    else if (!strcmp(EdgeWeightFormat, "LOWER_DIAG_COL"))
        WeightFormat = LOWER_DIAG_COL;
    else
        eprintf("Unknown EDGE_WEIGHT_FORMAT: %s", EdgeWeightFormat);
}

static void Read_EDGE_WEIGHT_SECTION()
{
    Node *Ni, *Nj;
    int i, j, W;

    CheckSpecificationPart();
    if (!FirstNode)
        CreateNodes();

    CostMatrix = (int *) calloc((size_t) Dimension * (Dimension - 1) / 2,
                                sizeof(int));
    Ni = FirstNode->Suc;
    do {
        Ni->C =
            &CostMatrix[(size_t) (Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
    }
    while ((Ni = Ni->Suc) != FirstNode);

    switch (WeightFormat) {
    case FULL_MATRIX:
        for (i = 1, Ni = FirstNode; i <= Dimension; i++, Ni = Ni->Suc) {
            for (j = 1; j <= Dimension; j++) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                if (j < i)
                    Ni->C[j] = W;
            }
        }
        break;
    case UPPER_ROW:
        for (i = 1, Ni = FirstNode; i < Dimension; i++, Ni = Ni->Suc) {
            for (j = i + 1, Nj = Ni->Suc; j <= Dimension;
                 j++, Nj = Nj->Suc) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                Nj->C[i] = W;
            }
        }
        break;
    case LOWER_ROW:
        for (i = 2, Ni = FirstNode->Suc; i <= Dimension; i++, Ni = Ni->Suc) {
            for (j = 1; j < i; j++) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                Ni->C[j] = W;
            }
        }
        break;
    case UPPER_DIAG_ROW:
        for (i = 1, Ni = FirstNode; i <= Dimension; i++, Ni = Ni->Suc) {
            for (j = i, Nj = Ni; j <= Dimension; j++, Nj = Nj->Suc) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                Nj->C[i] = W;
            }
        }
        break;
    case LOWER_DIAG_ROW:
        for (i = 1, Ni = FirstNode; i <= Dimension; i++, Ni = Ni->Suc) {
            for (j = 1; j <= i; j++) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                Ni->C[j] = W;
            }
        }
        break;
    case UPPER_COL:
        for (j = 2, Nj = FirstNode->Suc; j <= Dimension; j++, Nj = Nj->Suc) {
            for (i = 1; i < j; i++) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                Nj->C[i] = W;
            }
        }
        break;
    case LOWER_COL:
        for (j = 1, Nj = FirstNode; j < Dimension; j++, Nj = Nj->Suc) {
            for (i = j + 1, Ni = Nj->Suc; i <= Dimension;
                 i++, Ni = Ni->Suc) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                Ni->C[j] = W;
            }
        }
        break;
    case UPPER_DIAG_COL:
        for (j = 1, Nj = FirstNode; j <= Dimension; j++, Nj = Nj->Suc) {
            for (i = 1; i <= j; i++) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                Nj->C[i] = W;
            }
        }
        break;
    case LOWER_DIAG_COL:
        for (j = 1, Nj = FirstNode; j <= Dimension; j++, Nj = Nj->Suc) {
            for (i = j, Ni = Nj; i <= Dimension; i++, Ni = Ni->Suc) {
                if (!fscanint(ProblemFile, &W))
                    eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                Ni->C[j] = W;
            }
        }
        break;
    }
}

static void Read_FIXED_EDGES_SECTION()
{
    Node *Ni, *Nj, *N, *NPrev = 0, *NNext;
    int i, j, Count = 0;

    CheckSpecificationPart();
    if (!FirstNode)
        CreateNodes();
    if (!fscanint(ProblemFile, &i))
        i = -1;
    while (i != -1) {
        if (i <= 0 || i > Dimension)
            eprintf("FIXED_EDGES_SECTION: Node number out of range: %d",
                    i);
        fscanint(ProblemFile, &j);
        if (j <= 0 || j > Dimension)
            eprintf("FIXED_EDGES_SECTION: Node number out of range: %d",
                    j);
        if (i == j)
            eprintf("FIXED_EDGES_SECTION: Illegal edge: %d to %d", i, j);
        Ni = &NodeSet[i];
        Nj = &NodeSet[j];
        if (!FixEdge(Ni, Nj))
            eprintf("FIXED_EDGES_SECTION: Illegal fix: %d to %d", i, j);
        /* Cycle check */
        N = Ni;
        Count = 0;
        do {
            NNext = N->FixedTo1 != NPrev ? N->FixedTo1 : N->FixedTo2;
            NPrev = N;
            Count++;
        } while ((N = NNext) && N != Ni);
        if (N == Ni && Count != Dimension)
            eprintf("FIXED_EDGES_SECTION: Illegal fix: %d to %d", i, j);
        if (!fscanint(ProblemFile, &i))
            i = -1;
    }
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

    CheckSpecificationPart();
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
