#!/bin/sh

# this is to cut out a piece of the vstep and present this as a wiggle.


if [ $# -ne 5 ]
then
	echo "usage: vstepregion2wiggle.sh <vstepfile> <width> <chrom> <start> <end>"
	exit 1
fi


name=`basename $1`" - wiggle portion"

echo "track type=wiggle_0 name=\""$name"\" visibility=\"full\" description=\""$name"\" alwaysZero=\"on\""

#awk -v achr="$3" ' $1==achr { print $0 } ' $1

awk -v awidth="$2" -v achr="$3" -v astart="$4" -v aend="$5" '$1==achr && $2 >= astart && $2 <= aend { print $1"\t"$2"\t"$2+awidth"\t"$3 }' $1 

