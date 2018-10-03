#!/usr/bin/env python

import csv
import matplotlib.pyplot as plt
import numpy as np

infile = "/home/pzs/histone/HISTONE_DATA/uniq_exons_and_introns_NCBI35.txt"

reader = csv.reader(open(infile, "r"), delimiter="\t")

N = 400

canonical = [ 0 ] * N
alternative = [ 0 ] * N
intron = [ 0 ] * N

canonicallast = 0
alternativelast = 0
intronlast = 0

for row in reader:
	position = row[1][1:]
	# exon. Strip off first character to get canonical/alternative designation
	if row[2].startswith("E"):
		eitype = row[2][1:]
	else:
		eitype = "intron"
	if position == "last":
		if eitype == "canonical":
			canonicallast += 1
		elif eitype == "alternative":
			alternativelast += 1
		elif eitype == "intron":
			intronlast += 1
	else:
		if eitype == "canonical":
			filelist = canonical
		elif eitype == "alternative":
			filelist = alternative
		elif eitype == "intron":
			filelist = intron
		try:
			filelist[int(position)] += 1
		except IndexError:
			print position
			raise

canonical.append(canonicallast)
alternative.append(alternativelast)
intron.append(intronlast)


fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(N + 1)  # the x locations for the groups
width = 0.35       # the width of the bars

rects1 = ax.bar(ind+width, canonical, width)
rects2 = ax.bar(ind+width, alternative, width)
rects3 = ax.bar(ind+width, intron, width)

plt.show()
