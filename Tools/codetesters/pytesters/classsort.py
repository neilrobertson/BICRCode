#!/usr/bin/env python

class Sorter:
	def __init__(self, lookup):
		self.lookup = lookup
		
	def __call__(self, x, y):
		return cmp(self.lookup[x], self.lookup[y])

look = { 	"A" : 3,
			"B" : 5,
			"C" : 1,
			"D" : -1,
			"E" : 7 }
s = Sorter(look)

mylist = ["A","B","C","D","E"]
mylist.sort(s)

print mylist
