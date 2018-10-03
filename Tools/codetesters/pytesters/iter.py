#!/usr/bin/python

class HashWrapper:
	def __init__(self, inhash):
		self.hash = inhash
		
	def __iter__(self):
		return iter(self.hash)

	def __getitem__(self, item):
		return self.hash[item]

myhash = {	"A"	:	"hello",
			"B"	:	"goodbye",
			"C" :	"another" }
			
hw = HashWrapper(myhash)
for item in hw:
	hitem = hw[item]
	print item, hitem
