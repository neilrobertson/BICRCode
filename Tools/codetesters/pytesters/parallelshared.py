#!/usr/bin/env python

import pp
from string import lowercase

class MainList:
	def __init__(self):
		self.mainlist = []

	def pAppend(self, aitem):
		self.mainlist.append(aitem)
		print self.mainlist
	

ncpus = 2

job_server = pp.Server(ncpus, ppservers=())

ml = MainList()

jobs = [ job_server.submit(ml.pAppend ,(letter,), (), ()) for letter in lowercase ]

for job in jobs:
	job()
	
print ml.mainlist
