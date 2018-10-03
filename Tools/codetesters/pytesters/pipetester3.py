#!/usr/bin/env python

from subprocess import *
import csv

filename = "/home/pzs/workspace/scriptsandbox/largefile.txt"

cmd = "cat %s" % (filename)

p = Popen(cmd, shell=True, bufsize=100, stdout=PIPE)

reader = csv.reader(p.stdout, delimiter="\t")

for row in reader:
	chrom, coord, data = row
	if data > 0.8:
		print chrom, coord
