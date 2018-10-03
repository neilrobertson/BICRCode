#!/bin/bash

# a wget download script for pages that don't have a "download all" link

base="http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/data/"
extension=".bed.gz"
chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
types="RestingNucleosomes ActivatedNucleosomes"

for chrom in $chroms
do
	for dtype in $types
	do
		url=$base$dtype-chr$chrom$extension
		echo "wget $url"
		wget $url
	done
done
