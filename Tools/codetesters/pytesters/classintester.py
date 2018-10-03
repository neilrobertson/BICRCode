#!/usr/bin/env python

class top:
	def __init__(self, mylist):
		self.mylist = mylist
		
	def __contains__(self, item):
		return item in self.mylist
		
	def __len__(self):
		return len(self.mylist)

t = top([1,2,3,4,5])

print len(t)
