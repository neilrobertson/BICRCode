#!/usr/bin/env python

class top:
	def __init__(self):
		self.mylist = []

	def addItem(self, item):
		self.mylist.append(item)

	def __contains__(self, item):
		return item in self.mylist
		
class Inner:
	def __init__(self, left, right):
		self.left = left
		self.right = right
		
	def __eq__(self, other):
		return self.left == other.left and self.right == other.right

t = top()

item1 = Inner(3,4)

t.addItem(item1)

item2 = Inner(3,4)
item3 = Inner(4,4)

print item3 in t

