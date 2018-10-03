#!/bin/bash

# similar to cutall.sh

for i in *.txt
do
	newname=`basename $i .txt`-filtered.txt
	cut -f1,2,3 $i > processed/$newname
done
