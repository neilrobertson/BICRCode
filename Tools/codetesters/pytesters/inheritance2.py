#!/usr/bin/env python

class top:
	def __init__(self):
		self.name = "hello"
		self.tester()
		
class bottom(top):
	def __init__(self):
		top.__init__(self)
		print self.name

	def tester(self):
		print "inherited!"


b = bottom()

if b.__class__.__name__ == "top":
	print "yes"
