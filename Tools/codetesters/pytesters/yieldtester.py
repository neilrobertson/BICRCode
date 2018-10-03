#!/usr/bin/env python

def withoutYield(inlist):
	outlist = [ x ** 2 for x in inlist ]
	return outlist
	
def withYield(inlist):
	for x in inlist:
		yield x ** 2
		
		
startlist = range(10)
squares = withYield(startlist)
for square in squares:
	print square
