#!/usr/bin/env python

allgenes = "lhqsu80FNK.mart"
ourregions = "positions.txt"

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

regions = []
def getChromosomeRegions(chrfile):
	fh = open(chrfile, "r")
	for line in fh:
		line = line.replace(",", "")[:-1]
		chrom, rest = line.split(":")
		chrom = chrom[3:]
		start, end = rest.split("-")
		regions.append((chrom, start, end))
	return regions
	fh.close()
	
def getGeneRegions(genefile, regions):
	fh = open(genefile)
	header = fh.readline()
	for line in fh:
		line_sp = line[:-1].split(",")
		try:
			chrom = line_sp[0]
			start = line_sp[1]
			end = line_sp[2]
		except IndexError:
			breakpoint()
			continue
		if inRegions(chrom, start, end, regions):
			print line[:-1]
	fh.close()

def inRegions(chrom, start, end, regions):
	for region in regions:
		regchrom, regstart, regend = region
		if regchrom != chrom:
			continue
		if (start > regstart and start < regend) or (end > regstart and end < regend):
			return True
	return False
			

regions = getChromosomeRegions(ourregions)
getGeneRegions(allgenes, regions)

