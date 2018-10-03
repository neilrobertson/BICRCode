#!/usr/bin/env python

from breakpoint import *

testdict = 	{	"A"	:	{	"1"	: [1],
							"2" : [2]	},
				"B"	:	{	"1"	: [3],
							"2" : [4]	},
				"C"	:	{	"1"	: [5],
							"2"	: [6],	}}

							
def rotate(indict):
	breakpoint()
	rotated = {}
	for key1 in indict:
		innerdict = indict[key1]
		for key2 in innerdict:
			data = innerdict[key2]
			rotatedinner = rotated.setdefault(key2, {})
			rotatedinner[key1] = data
	return rotated
	
print rotate(testdict)
