#!/usr/bin/env python

class Wrapper(dict):
	def __init__(self):
		pass
		
w = Wrapper()
w["hello"] = 3
print w["hello"]
