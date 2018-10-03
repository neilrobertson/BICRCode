#!/usr/bin/env python

class ValList(list):
	def __init__(self):
		pass
		
	def type(self):
		return "vallist"

	def __sub__(self, other):
		newvl = ValList()
		if other.type() == "val":
			for i in range(len(self)):
				thispoint = self[i]
				newval = thispoint - other
				newvl.append(newval)
		elif other.type() == "vallist":
			for i in range(len(self)):
				thispoint = self[i]
				otherpoint = other[i]
				newval = thispoint - otherpoint
				newvl.append(newval)
		else:
			raise TypeError
		return newvl

class Val:
	def __init__(self, val):
		self.val = val
		
	def __sub__(self, other):
		if other.type() != "val":
			raise TypeError, "can only deduct a Val from another Val"
		newval = self.val - other.val
		newt = Val(newval)
		return newt
		
	def type(self):
		return "val"

	def __repr__(self):
		return "Val(%s)" % (self.val)
		
v1 = Val(1)
v2 = Val(2)
v3 = Val(3)
v4 = Val(4)

vl1 = ValList()
vl1.extend([ v4, v3 ])

vl2 = ValList()
vl2.extend([ v1, v2 ])


if isinstance(v1, Val):
	print "yes!"
