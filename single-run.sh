#!/bin/sh
#
# This script takes a timeout in seconds and an input file
# and lets the release version of irp run on the input file;
# output is written to results/<name based on commit or tag>/out_*.txt
#
# Make sure you have compiled with
#    cargo build --release
#
# We could, for example do
#    cargo build --release && echo instances/*/*/*.dat | xargs -P 4 -n 1 ./single-run.sh 5

timeout="$1"
test -n "$timeout" || exit 1
filename="$2"
test -n "$filename" || exit 2

resultdir="results/$(git describe --tags --always)"
mkdir -p "$resultdir"

basename="$(basename "$filename")"
logfile="$resultdir/$(echo "$basename" | sed -e 's/^\(.*\)\.dat$/out_\1.log/')"
resultfile="$resultdir/$(echo "$basename" | sed -e 's/^\(.*\)\.dat$/out_\1.txt/')"

echo "Running $timeout seconds for input $filename and writing to $resultfile"
timeout $timeout ./target/release/plaisir "$filename" 2>"$resultfile" >"$logfile"
echo "Finished $filename"
