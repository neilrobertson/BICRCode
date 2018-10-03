#!/usr/bin/python

def coordcmp(x1, x2):
	return cmp(x1[1], x2[1])

mylist = [('d',5,5), ('a',2,3), ('b',4,2), ('c',1,3) ]

mylist.sort(coordcmp)

print mylist
