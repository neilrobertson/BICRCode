#!/usr/bin/env python

def mycmp(a, b):
	return cmp(a, b)

mylist = [6,1,61,7,1,34,16]
mylist.sort(mycmp)

print mylist
