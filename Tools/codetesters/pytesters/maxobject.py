#!/usr/bin/env python

class ValueObject:
	def __init__(self, value):
		self.value = value
		
	def __ge__(self, other):
		return self.value > other.value
		
	def __repr__(self):
		return "ValueObject(%s)" % (self.value)

		
o1 = ValueObject(3.0)
o2 = ValueObject(3.6)
o3 = ValueObject(6.0)

mylist = [o1,o2,o3]

print max(mylist)
