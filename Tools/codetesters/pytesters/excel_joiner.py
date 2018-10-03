#!/usr/bin/env python

import csv

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

incsv = "sreenivas.csv"

def median(vallist):
	num_vals = len(vallist)
	vallist.sort()
	if num_vals % 2 == 1: # odd
		index = (num_vals - 1) / 2
		return vallist[index]
	else: # even
		index = num_vals / 2
		return (vallist[index] + vallist[index - 1]) / 2

def oneGroup(line_sp, startindex, alldata, gname):
	name = line_sp[startindex]
	id = line_sp[startindex + 1].strip('"')
	if not name or not id:
		return
	unid = "%s--%s" % (name, id)
	data = line_sp[startindex + 2]
	this_map = alldata.setdefault(unid, {})
	this_map[gname] = data

alldata = {}

fh = open(incsv, "r")
header = fh.readline()
for line in fh:
	line_sp = line[:-1].split(",")
	oneGroup(line_sp, 0, alldata, "group1")
	oneGroup(line_sp, 4, alldata, "group2")
	oneGroup(line_sp, 8, alldata, "group3")
	oneGroup(line_sp, 12, alldata, "group4")

index = alldata.keys()
index.sort()

lines = []
line1 = [ "name", "id", "data1", "data2", "data3", "data4", "median" ]
lines.append(line1)
for ind in index:
	data = alldata[ind]
	id, name = ind.split("--")
	points = [data.get("group1", ""), data.get("group2", ""), data.get("group3", ""), data.get("group4", "")]
	intpoints = []
	for point in points:
		try:
			intpoints.append(float(point))
		except:
			pass
	if not intpoints:
		breakpoint()
	av_val = median(intpoints)
	line = [ id, name ] + points + [ str(av_val) ]
	lines.append(line)
	

writer = csv.writer(open("output.csv", "wb"))
writer.writerows(lines)
