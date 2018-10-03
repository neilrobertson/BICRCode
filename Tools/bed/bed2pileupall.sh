if [ $# -ne 3 ]
then
	echo "usage: bed2pileupall.sh <beddirectory> <extension> <outdir>"
	exit 1
fi

b2pcmd="python /home/pzs/repository/shared/tools/bed/bed2pileupwiggle.py"


for i in $1/*.$2
do
	cmd="$b2pcmd $i $3"
	echo $cmd
	$cmd
done
