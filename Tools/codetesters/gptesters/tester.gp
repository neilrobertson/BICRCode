set terminal png
set output 'tester.png'
plot 	\
		'tester.data' using 1:2 title "y" with lines,\
		'tester.data' using 1:3 title "z" with lines
