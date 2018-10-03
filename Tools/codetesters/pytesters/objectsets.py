#!/usr/bin/env python

from sets import Set

class Tester:
	def __init__(self, i1, i2):
		self.i1 = i1
		self.i2 = i2

s1 = Set()

t1 = Tester(1,2)
t2 = Tester(1,2)

s1.add(t1)
s1.add(t2)

print s1
