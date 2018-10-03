#!/bin/bash

# similar to filterall.sh

for i in *binned.csv
do
	newfile=`basename $i .csv`2.csv
	cut -f1,2,4 $i > $newfile
	mv $newfile $i
done
