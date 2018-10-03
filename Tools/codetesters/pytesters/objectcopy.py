#!/usr/bin/env python

from copy import deepcopy

class TesterObject:
	def __init__(self):
		self.mylist = [1]
		
	def appendItem(self, item):
		self.mylist.append(item)


to = TesterObject()
to2 = deepcopy(to)

to2.appendItem(3)

print to.mylist
print to2.mylist
