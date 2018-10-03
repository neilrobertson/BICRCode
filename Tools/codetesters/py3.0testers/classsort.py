#!/usr/bin/env python3.0

class Obj:
	def __init__(self, val):
		self.val = val
		
	def __lt__(self, other):
		return self.val < other.val

mylist = []

o1 = Obj(3)
o2 = Obj(4)
o3 = Obj(10)

mylist = [ o2, o3, o1 ]
 
mylist.sort()
