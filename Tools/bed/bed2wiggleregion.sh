#!/bin/sh

# this is to cut out a piece of the bed and present this as a wiggle.


if [ $# -ne 6 ]
then
	echo "usage: bed2wiggleregion.sh <bedfile> <chrom> <start> <end> <headername> <outdir>"
	exit 1
fi


name=`basename $1`" - $5"

outfile=$6/$1.$5

# echo "track type=wiggle_0 name=\""$name"\" visibility=\"full\" description=\""$name"\" alwaysZero=\"on\" \r" > $outfile
# echo "browser position $2:$3-$4\r" >> $outfile

#awk -v achr="$3" ' $1==achr { print $0 } ' $1

awk -v achr="$2" -v astart="$3" -v aend="$4" '$1==achr && $2 <= aend && $3 >= astart { print }' $1 >> $outfile

