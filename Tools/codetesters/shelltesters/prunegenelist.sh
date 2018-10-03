#!/bin/bash

if [ $# -ne 1 ]
then
	echo "usage: prunegenelist.sh <genelist>"
	exit 1
fi

cut -f 2,3,4,5,6,7,8 genes-human-CD4-NCBIm35i.csv | sort -k1 -k4 -k7 | uniq
