#!/usr/bin/env python

class CoordLookup:
	def __init__(self):
		self._coords = {}
		self._names = {}
		self._strands = {}
		
	def setFromFile(self, fname, fsettings):
		print "setting up coordinate lookup for file", fname
		import csv
		idindex = fsettings["id"]
		if "name" in fsettings:
			nameindex = fsettings["name"]
		chromindex = fsettings["chrom"]
		if "strand" in fsettings:
			strandindex = fsettings["strand"]
		startindex = fsettings["start"]
		endindex = fsettings["end"]
		reader = csv.reader(open(fname), delimiter="\t")
		for row in reader:
			try:
				id_ = row[idindex]
				if "name" in fsettings:
					name = row[nameindex]
				else:
					name = ""
				if "strand" in fsettings:
					strand = row[strandindex]
				chrom = row[chromindex]
				start = int(row[startindex])
				end = int(row[endindex])
			except IndexError:
				print row
				raise
			if "strand" in fsettings:
				self._strands[id_] = strand 
			if not chrom.startswith("chr"):
				chrom = "chr" + chrom
			self._coords[id_] = (chrom, start, end)
			if name != "":
				self._names[id_] = name

	def getStrand(self, gid):
		return self._strands[gid]

	def getCoord(self, gid):
		return self._coords[gid]

	def getName(self, gid):
		return self._names[gid]

	def __iter__(self):
		return iter(self._coords)
