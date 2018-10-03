#!/bin/sh

chroms="1 
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
X
Y"

basefile="/home/pzs/codetesters/perltesters/gtacwiggle2.txt"
outdir="/home/pzs/codetesters/shelltesters/output/"

for i in $chroms
do
	echo "working at chrom" $i
	outfile=$outdir"chrom"$i".txt"
	echo 'browser hide all' > $outfile
	echo 'browser full altGraph' >> $outfile
	echo 'track type=wiggle_0 name="GTAC: chr'$i'" description="GTAC: chr'$i'" visibility=full color=0,0,0 viewLimits=0:2 yLineMark=0 yLineOnOff=on useScore=0 priority=0 offset=0 autoscale=on gridDefault=off maxheightPixels=128:128:5 graphType=bar' >> $outfile
	grep "^chr$i[^0-9]" $basefile >> $outfile
done
