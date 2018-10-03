#!/bin/bash

# runs the UCSC tool to turn a wig into a bigwig, but also generates a standard header file

baseurl="http://cs-analysis.clinmed.gla.ac.uk/bigwigs"
localbase="/var/www/bigwigs"
w2bw="$HOME/bin/wigToBigWig"

if [ $# -ne 3 ]
then
	echo "usage: wig2BigWig.sh <wigfile> <buildname> <outname>"
	exit 1
fi

chrmfile="$HOME/mount/publicdata/$2/chrmSizes.$2"

# the filebase is everything before the first "."
filebase=`basename $1 | cut -d"." -f1`
samplename=`basename $filebase`

echo "filebase " $filebase
echo "samplename " $samplename

outdir=$localbase/$3
echo "***dumping into directory" $outdir
if ! [ -a $outdir ]
then
	mkdir $outdir
fi
outbw=$outdir/$filebase.bw
outheader=$outbw.txt

if [ -a $outbw ]
then
	echo "***file already exists: $outbw"
	exit 0
fi

cmd="$w2bw $1 $chrmfile $outbw"
echo $cmd
$cmd

if [ $? -ne 0 ]
then
	exit 1
fi

url=$baseurl/$3/$filebase.bw

echo "track type=bigWig name=\"$samplename\" description=\"$samplename\" bigDataUrl=$url" > $outheader
