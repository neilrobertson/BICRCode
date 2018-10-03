#!/usr/bin/env python

threshold = 1000

class Tester:
	def __init__(self, start, end):
		self.start = start
		self.end = end
		self.mid = self.end - self.start
		
	def __eq__(self, other):
		if self.mid == other.mid:
			return True
		if abs(self.mid - other.mid) < threshold:
			return doesOverlap(self.start, self.end, other.start, other.end)
		return False
			
def doesOverlap(start, end, lb, ub):
	# outside - not in bin at all
	if ub < start or lb > end:
		return False
	# completely inside bin
	if ub <= end and lb >= start:
		return True
	if ub >= end and lb <= start:
		return True
	# to left of bin
	if lb < start and ub > start and ub < end:
		return True
	# to right of bin
	if ub > end and lb < end and lb > start:
		return True
			
t1 = Tester(100, 200)

mylist = [t1]



t2 = Tester(201, 300)

print t1 == t2
