#!/bin/sh

# Configuration
TIMEOUT=300
INSTANCES='L_abs1n50_4_L L_abs5n200_4_H L_abs3n50_4_L L_abs10n200_5_H'

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
