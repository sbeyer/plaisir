#!/bin/sh

# Configuration
VERIFIER="../dimacs-irp-verifier/verify.py"

exitcode=0
for file in "$@"
do
	# Instance name
	i="$(basename "$file" .txt | sed 's/^out_//')"

	# Input file
	input="$(git ls-files "instances/**$i*.dat")"

	# Log file for verification
	verifylog="/tmp/$i.verify.log"

	"$VERIFIER" "$input" "$(dirname "$file")" >$verifylog
	if test "$?" -ne 0
	then
		echo
		echo "Verification of $i failed, see $verifylog"
		tail -n 2 "$verifylog" | sed 's/^/>>> /'
		exitcode=1
	fi
done

exit $exitcode
