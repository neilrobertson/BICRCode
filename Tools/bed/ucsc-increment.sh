#!/bin/sh

if [ $# -ne 1 ]
then
	echo "usage: ucsc-increment.sh <bedfile>"
	exit 1
fi

awk -F"\t" '{ if (/^#/)
       print
     else
       { 
         printf $1 "\t"  $2+1 "\t"   $3+1 "\t";
	 for (i=4;i<=NF;i++)
         {
           printf("%s\t",$i)
         }
	 printf("\n")
        };
      }' $1

