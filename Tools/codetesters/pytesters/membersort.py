#!/usr/bin/env python

from operator import attrgetter

class Tester:
	def __init__(self, innum):
		self._num = innum
		
	def __repr__(self):
		return str(self._num)
		
tlist = [ Tester(3), Tester(4), Tester(1) ]

tlist.sort(key=attrgetter("_num"))

print tlist
