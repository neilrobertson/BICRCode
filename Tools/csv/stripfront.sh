for i in *.txt
mv $i `echo $i | cut -d"_" -f 2`
