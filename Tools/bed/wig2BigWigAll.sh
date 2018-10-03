#!/bin/bash

if [ $# -ne 4 ]
then
	echo "usage: wig2BigWigAll.sh <dir> <extension> <buildname> <outname>"
	exit 1
fi

w2bws=$HOME/repository/shared/tools/bed/wig2BigWig.sh

for i in $1/*$2
do
	cmd="$w2bws $i $3 $4"
	echo $cmd
	$cmd

	if [ $? -ne 0 ] 
	then
		echo "problem at file: $i"
		exit 1
	fi
done

