#!/bin/sh

# Configuration
TIMEOUT=120
INSTANCES='S_abs4n25_2_L3 S_abs5n10_2_H6 S_abs2n35_4_L3 S_abs4n15_5_H6'
VERIFIER="../dimacs-irp-verifier/verify.py"

resultdir="results/$(git describe --tags --always)"

# Build
cargo build --release || exit

# Find instances
FILES=""
for i in $INSTANCES
do
	FILES="$FILES $(git ls-files "instances/**$i*.dat")"
done

# Solve
echo "$FILES" | xargs -P 4 -n 1 ./single-run.sh "$TIMEOUT"

# Verify
for file in $FILES
do
	i="$(basename "$file" .dat)"
	verifylog="/tmp/$i.verify.log"
	"$VERIFIER" "$file" "$resultdir/" >$verifylog ||
		echo "Verification of $i failed, see $verifylog"
done

# Obtain results
./scrape-results.py "$resultdir"/*.txt
