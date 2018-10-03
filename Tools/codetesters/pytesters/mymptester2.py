#!/usr/bin/env python

from __future__ import print_function
from multiprocessing import Pool, current_process
import random
import time
from string import lowercase

ncpus = 3

def pipeline1(instr):
	time.sleep(random.random())
	return instr + "-p1-" + current_process().name

def pipeline2(instr):
	time.sleep(random.random())
	return instr + "-p2-" + current_process().name

pool = Pool(processes=ncpus)

before = lowercase
afterp1 = pool.map(pipeline1, before)
afterp2 = pool.map(pipeline2, afterp1)

map(print, afterp2)
