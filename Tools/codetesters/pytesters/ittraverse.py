#!/usr/bin/env python

from cintervaltree import *

class IntervalCount:
	def __init__(self):
		self._count = 0
		
	def __call__(self, iv):
		if iv.get_interval().value >= 2:
			self._count += 1
			
	def getValue(self):
		return self._count

def printValue(iv):
	print iv.get_interval().value

def setZero(iv):
	iv.get_interval().value = 0

it = IntervalTree()
i1 = Interval(10, 20, value=0)
i2 = Interval(30, 40, value=3)
i3 = Interval(40, 50, value=3)
i4 = Interval(50, 60, value=3)
it.insert_interval(i1)
it.insert_interval(i2)
it.insert_interval(i3)
it.insert_interval(i4)


newcount = IntervalCount()
it.traverse(newcount)
print newcount.getValue()
it.traverse(setZero)

newcount = IntervalCount()
it.traverse(newcount)
print newcount.getValue()
