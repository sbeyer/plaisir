#!/bin/sh

# Configuration
TIMEOUT=120
INSTANCES='S_abs4n25_2_L3 S_abs5n10_2_H6 S_abs2n35_4_L3 S_abs4n15_5_H6'

resultdir="results/$(git describe --tags --always)"

# Build
cargo build --release || exit

# Find instances
INPUT_FILES=""
for i in $INSTANCES
do
	INPUT_FILES="$INPUT_FILES $(git ls-files "instances/**$i*.dat")"
done

# Solve
echo "$INPUT_FILES" | xargs -P 4 -n 1 ./single-run.sh "$TIMEOUT"

# Collect output files
OUTPUT_FILES=""
for i in $INSTANCES
do
	OUTPUT_FILES="$OUTPUT_FILES $resultdir/out_$i.txt"
done

# Verify
./verify.sh $OUTPUT_FILES || exit

# Obtain results
./scrape-results.py $OUTPUT_FILES
