#!/usr/bin/env python

class SampleException(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


# apparently this is done differently in python 2.6 (and therefore, I guess, python 3.0)
try:
	raise SampleException("some text")
except SampleException, e:
	print e.value
