#!/bin/sh

# Configuration
VERIFIER="../dimacs-irp-verifier/verify.py"

resultdir="results/$(git describe --tags --always)"
test -n "$1" && resultdir="$1"
test -d "$resultdir" || ( echo "Result directory $resultdir does not exist" ; exit 1 ) || exit

# Find instances
FILES=""
for file in "$resultdir"/*.txt
do
	i="$(basename "$file" .txt | sed 's/^out_//')"
	FILES="$FILES $(git ls-files "instances/**$i*.dat")"
done

# Verify
exitcode=0
for file in $FILES
do
	i="$(basename "$file" .dat)"
	verifylog="/tmp/$i.verify.log"
	"$VERIFIER" "$file" "$resultdir/" >$verifylog
	if test "$?" -ne 0
	then
		echo "Verification of $i failed, see $verifylog"
		exitcode=1
	fi
done

exit $exitcode
