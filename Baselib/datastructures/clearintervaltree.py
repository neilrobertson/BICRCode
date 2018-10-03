#!/usr/bin/env python

# object to be used with CIntervalTree's traverse method
# sets all interval nodes to have the specified value
class ClearIntervalTree:
	def __init__(self, setval):
		self._setval = setval
		
	def __call__(self, iv):
		iv.get_interval().value = self._setval
