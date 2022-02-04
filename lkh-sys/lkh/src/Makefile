TREE_TYPE = TWO_LEVEL_TREE
# TREE_TYPE = THREE_LEVEL_TREE
# TREE_TYPE = ONE_LEVEL_TREE

# CC = gcc
IDIR = include
ODIR = obj
CFLAGS = -O3 -flto -Wall -I$(IDIR) -D$(TREE_TYPE) -g -DNDEBUG

_DEPS = Delaunay.h GainType.h GeoConversion.h Hashing.h      \
        Heap.h libLKH.h LKH.h Segment.h Sequence.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = Activate.o AddCandidate.o AddExtraCandidates.o                  \
       AddTourCandidates.o AdjustCandidateSet.o AdjustClusters.o       \
       AllocateStructures.o Ascent.o                                   \
       Best2OptMove.o Best3OptMove.o Best4OptMove.o Best5OptMove.o     \
       BestKOptMove.o Between.o Between_SL.o Between_SSL.o             \
       BridgeGain.o BuildKDTree.o C.o CandidateReport.o                \
       ChooseInitialTour.o Connect.o CreateCandidateSet.o              \
       CreateDelaunayCandidateSet.o CreateNNCandidateSet.o             \
	   CreateQuadrantCandidateSet.o     \
       Delaunay.o Distance.o eprintf.o ERXT.o       \
       Excludable.o Exclude.o FindTour.o FixedOrCommonCandidates.o     \
       Flip.o Flip_SL.o Flip_SSL.o FreeStructures.o        \
       fscanint.o Gain23.o GenerateCandidates.o                        \
       GeoConversion.o GetTime.o GreedyTour.o Hashing.o Heap.o         \
       IsCandidate.o \
       IsPossibleCandidate.o KSwapKick.o LinKernighan.o LKHmain.o      \
       Make2OptMove.o Make3OptMove.o Make4OptMove.o Make5OptMove.o     \
       MakeKOptMove.o MergeTourWithBestTour.o MergeWithTourIPT.o       \
       Minimum1TreeCost.o MinimumSpanningTree.o NormalizeNodeList.o    \
       NormalizeSegmentList.o OrderCandidateSet.o PatchCycles.o        \
       printff.o \
       Random.o \
	   RecordBestTour.o \
	   RecordBetterTour.o RemoveFirstActive.o                          \
       ResetCandidateSet.o RestoreTour.o SegmentSize.o Sequence.o      \
       SFCTour.o \
       Statistics.o StoreTour.o SymmetrizeCandidateSet.o               \
       TrimCandidateSet.o \
       LKH.o

all:
	$(MAKE) LKH

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR):
	mkdir -p "$(ODIR)"

.PHONY:
	all clean

LKH: $(ODIR) $(OBJ) $(DEPS)
	ar rcs ../libLKH.a $(OBJ)

clean:
	rm -rf $(ODIR) ../libLKH.a *~ ._* $(IDIR)/*~ $(IDIR)/._*

