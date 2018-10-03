#!/bin/bash

# this also removes chrM and trims to the first three columns
# this isn't quite satisfactory, because chromosome 2 comes after chromosome 10. 
# however, in this case, what's most important is that each chromosome is separated
# if we use a numerical sort for chromosomes, it mixes X and Y together.

if [ $# -ne 1 ]
then
	echo "usage: sortbed.sh <beddirectory>"
	exit 1
fi

for i in $1/*.bed
do
	echo "working at file:" $i
	sort -k 1.4,1.5 -k 2,2n $i | grep -v chrM | cut -f1,2,3 > $i.sort
done
