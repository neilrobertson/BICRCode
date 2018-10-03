#!/bin/sh

# this is to cut out a piece of the vstep and present this as a wiggle.


if [ $# -ne 2 ]
then
	echo "usage: vstep2wiggle.sh <vstepfile> <width>"
	exit 1
fi


name=`basename $1`" - wiggle"

echo "track type=wiggle_0 name=\""$name"\" visibility=\"full\" description=\""$name"\" alwaysZero=\"on\""

#awk -v achr="$3" ' $1==achr { print $0 } ' $1

awk -v awidth="$2" '{ print $1"\t"$2"\t"$2+awidth"\t"$3 }' $1 

