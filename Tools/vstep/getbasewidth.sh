#!/bin/sh

for i in ../zhao09HDAC/*vstep.txt
do
	base=`basename $i .txt`
	echo $base `grep variableStep $i | head -1 | cut -d' ' -f 3 | cut -d= -f 2`
done
