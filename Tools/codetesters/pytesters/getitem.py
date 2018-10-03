#!/usr/bin/env python

class ListWrapper:
	def __init__(self, inlist):
		self.mylist = inlist
		
	def __getitem__(self, index):
		return self.mylist[index]
		
testlist = [1,2,3,4]

lw = ListWrapper(testlist)

print lw[2]
