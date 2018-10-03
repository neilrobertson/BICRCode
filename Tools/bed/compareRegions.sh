#!/bin/dash

for i in $1/*
do
    for j in $1/*
    do
        target=${i#"./"}
        control=${j#"./"}

        if [ $target = "vs" ]
        then
            continue
        fi
	if [ $target = "compare.sh" ]
        then
            continue
        fi
        if [ $control = "vs" ]
        then
            continue
        fi
        if [ $control = "compare.sh" ]
        then
            continue
        fi
        if [ $target = $control ]
        then
            continue
        fi

        echo "****"
        echo $target" vs "$control
        echo "****";

	python $HOME/mount/repository/shared/baselib/bed/overlap.py -b $target/regions.csv -a $control/regions.csv -p 100
    done
done


