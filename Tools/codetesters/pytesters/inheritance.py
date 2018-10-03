#!/usr/bin/env python

class top:
	def __init__(self):
		self.name = "hello"
		
class bottom(top):
	def __init__(self):
		top.__init__(self)
		print self.name


b = bottom()

if isinstance(b, top):
	print "yup"
