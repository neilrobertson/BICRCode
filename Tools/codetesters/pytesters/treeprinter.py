#!/usr/bin/env python


def printTree(tree,depth):
	if ( tree == None or len(tree.keys()) == 0 ):
		for i in range(depth):
			print "\t",
		print "-"
	else:
		for key in tree.keys():
			for i in range(depth):
					print "\t",
			print key
			printTree(tree[key],depth+1)

tester = {	"a"	:	{	"aa"	:	1,
						"bb"	:	2	},
			"b" : 	{	"cc"	:	3,
						"dd"	:	4,
						"ee"	:	5,	},
			"c"	:	{	"ff"	:	6	}
		}
printTree(tester,0)
						
