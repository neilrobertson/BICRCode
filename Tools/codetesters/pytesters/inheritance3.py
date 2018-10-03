#!/usr/bin/env python

class MyList(list):
	def __init__(self):
		pass
		
	def setSelf(self, inlist):
		self = inlist
		print self
		
	def first(self):
		return self[0]


ml = MyList()
print dir(ml)
