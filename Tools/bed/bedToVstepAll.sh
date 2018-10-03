#!/bin/bash

# strips off <bedextension> and adds <vstepwidth>.vstep as the output extension

b2v="$HOME/repository/shared/baselib/bed/bedToVstep.py"

if [ $# -ne 4 ]
then
	echo "usage: bedToVstepAll.sh <bedextension> <beddirectory> <outdir> <vstepwidth>"
	exit 1
fi

for i in $2/*$1
do
	cd $2
	base=`basename $i`
	out=$3/$base.$4.vstep
	echo python $b2v -i $i -o $out -s $4
	python $b2v -i $i -o $out -s $4
done
