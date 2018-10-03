#!/usr/bin/env python

import pp
import random
import time
from string import lowercase

ncpus = 3

class Pipeline:
	def __init__(self):
		pass
	
	def pp1extra(self, jdetail):
		return jdetail + " phase 1 complete"

	def pipeline1(self, waittime, jdetail):
		bothstages = "abcfgijkmnopqrsuvwz"
		mycore = open("/proc/%i/stat" % os.getpid()).read().split()[38]
		print "I'm doing stuff!", "--core:", mycore, "--job detail:", jdetail
		jfinished = self.pp1extra(jdetail)
		if jdetail not in bothstages:
			jfinished += " (done)"
	#	time.sleep(waittime)
		return jfinished

	def pipeline2(self, waittime, jdetail):
		mycore = open("/proc/%i/stat" % os.getpid()).read().split()[38]
		print "I'm doing stuff!", mycore, jdetail
		jfinished = jdetail + " phase 2 complete"
	#	time.sleep(waittime)
		return jfinished

	

job_server = pp.Server(ncpus, ppservers=())

ofh = open("output.txt", "w")

po = Pipeline()
# outer loop. The equivalent of stepping through the input file
batch1 = []
batch2 = []
batchsize = 10
for letter in lowercase:
	batch1.append(letter)
	if len(batch1) >= batchsize:
		# pass this batch for processing
		jobs1 = [ job_server.submit(po.pipeline1 ,(random.random(), jinput), (), ("os", "time",)) for jinput in batch1 ]
		for job in jobs1:
			if not job().endswith("(done)"):
				batch2.append(job())
		batch1 = []
	if len(batch2) >= batchsize:
		jobs2 = [ job_server.submit(po.pipeline2 ,(random.random(), jinput), (), ("os", "time",)) for jinput in batch2 ]
		for job in jobs2:
			print >> ofh, job()
		batch2 = []
			
