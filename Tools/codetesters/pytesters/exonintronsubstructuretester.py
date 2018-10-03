#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import pylab as plt


exons = plt.array([ [1,1,2], [3,2,3], [4,2,2], [4,5,1] ])
introns = plt.array([ [0.1, 0.1, 0.2, 0.5], [0.1, 0.4, 0.2, 0.4], [0.4, 0.1, 0.5, 0.2] ])

# width of each exon plot section 
width = 0.5

 
# exonxs = []
# intronxs = []
# exonlength = len(exons[0])
# intronlength = len(introns[0])
# exonbinlength = width / (exonlength - 1)
# intronbinlength = width / (intronlength - 1)
# currentexonx = 0
# currentintronx = width
# for i in range(len(exons)):
# 	for j in range(exonlength):
# 		exonxs.append(currentexonx)
# 		currentexonx += exonbinlength
# 	currentexonx -= exonbinlength
# 	currentexonx += width
# 	for j in range(intronlength):
# 		intronxs.append(currentintronx)
# 		currentintronx += intronbinlength
# 	currentintronx -= intronbinlength
# 	currentintronx += width
# 
# print exonxs
# print intronxs
# 
# flatexons = exons.flatten()
# flatintrons = introns.flatten()
# 
# plt.plot(exonxs, flatexons)
# plt.plot(intronxs, flatintrons)

exonlength = len(exons[0])
intronlength = len(introns[0])

for i in range(len(exons)):
	exon = exons[i]
	if i <= len(introns) - 1:
		plotintrons = True
	else:
		plotintrons = False
	exonstart = 2 * i * width 
	exonend = (2 * i * width) + width
	exonxs = plt.linspace(exonstart, exonend, exonlength)
	exlines = plt.plot(exonxs, exon, "blue")
	if plotintrons:
		intron = introns[i]
		intronstart = (2 * i * width) + (width)
		intronend = (2 * i * width) + (2 * width)
		intronxs = plt.linspace(intronstart, intronend, intronlength)
		inlines = plt.plot(intronxs, intron, "red")

def makeExonIntronXAxisLabels(num):
	labels = []
	# we only go from 1 to num because we're going to add two more for the last entries
	for i in range(1, num):
		elabel = "e" + str(i)
		ilabel = "i" + str(i)
		labels.extend([elabel, ilabel])
	labels = labels[:-1]
	labels.extend([ "Ilast", "Elast"])
	return labels

xlabels = makeExonIntronXAxisLabels(len(exons))
print len(xlabels) * width, width
xtickpositions = plt.linspace(0, (len(xlabels)-1)*width, len(xlabels)) + width / 2 
print xlabels
print xtickpositions
plt.xticks(xtickpositions,  xlabels, rotation="45")
plt.legend((exlines, inlines), ("Exons", "Introns"))
plt.savefig("substructure.png")



def plotExonIntronSubstructure(basedir, eiscomposite):
	plotpaths = {}
	for treatment in eiscomposite.treatments:
		# each value dictionary is indexed by the value in subplots
		subplots, exonvalues, intronvalues = eiscomposites.getValues()
		width = 0.5
		f = pylab.figure()
		pylab.subplots_adjust(hspace=0.7)
		name = "%s-%s" % (eiscomposite.name, treatment)
		filename = "%s/%s.png" % (basedir, name)
		print "dumping exon/intron partition plot:", filename
		numplots = len(subplots)
		thisplotnum = 1
		for spname in subplots:
			s = pylab.subplot(numplots, 1, thisplotnum)
			thisplotnum += 1
			exons = exonvalues[spname]
			introns = intronvalues[spname]
			exonwidth = len(exons[0])
			intronwidth = len(introns[0])
			for i in range(len(exons)):
				exon = exons[i]
				# number of introns is allowed to be smaller
				if i <= len(introns) - 1:
					plotintrons = True
				else:
					plotintrons = False
				assert len(exon) == exonlength
				exonstart = 2 * i * width 
				exonend = (2 * i * width) + width
				exonxs = plt.linspace(exonstart, exonend, exonlength)
				exlines = plt.plot(exonxs, exon, "blue")
				if plotintrons:
					intron = introns[i]
					intronstart = (2 * i * width) + (width)
					intronend = (2 * i * width) + (2 * width)
					intronxs = plt.linspace(intronstart, intronend, intronlength)
					inlines = plt.plot(intronxs, intron, "red")
			xlabels = makeExonIntronXAxisLabels(len(exons))
			numparts = len(labels)
			xtickpositions = plt.linspace(0, (numparts-1)*width, numparts) + width / 2 
			pylab.xticks(xtickpositions, xlabels, rotation="45")

class ExonIntronSubstructureComposite:
	def __init__(self, name):
		self.name = name
		# this is substructure composite, so these are lists of lists
		# each sublist is the substructure of one exon/intron
		self.exonbins = []
		self.intronbins = []
		self.lastexonbins = None
		self.lastintronbins = None
		self.numbins = 0
		self.exongran = 0
		self.introngran = 0
		
	# numbins is the maximum number of bins we'll need in each category
	# exongran and introngran are the granularity of the substructure
	# - how many bins to each exon and intron
	def initBins(self, numbins, exongran, introngran):
		for i in range(numbins):
			exonsubstructure = []
			intronsubstructure = []
			for j in range(exongran):
				subebin = AbstractSpreadBin(j, j+1)
				subebin.setOwner(self)
				exonsubstructure.append(subebin)
			for j in range(introngran):
				subibin = AbstractSpreadBin(i, i+1)
				subibin.setOwner(self)
				intronsubstructure.append(subibin)
			self.exonbins.append(ebin)
			self.intronbins.append(ibin)
		self.lastexonbin = AbstractSpreadBin(999, 1000)
		self.lastintronbin = AbstractSpreadBin(999, 1000)
		self.lastexonbin.setOwner(self)
		self.lastintronbin.setOwner(self)
		self.numbins = numbins
			
	# subbinlist is the composites list of bins - we're merging into here
	# the mergelist is what we're merging. 
	# strand denotes whether we need to reverse the list first, to take account of strand
	def mergeSubBins(self, subbinlist, mergelist, strand):
		assert len(subbinlist) == len(mergelist)
		if strand == -1:
			mergelist = mergelist[::-1]
		for i in range(len(mergelist)):
			subbin = subbinlist[i]
			mergebin = mergelist[i]
			subbin.mergeBin(mergebin)

	def addGene(self, gene):
		# a list of pairs
		# each pair is an exon/intron label and a list of the sub-bins
		fbins = gene.exonintronsubbins
		eindex = 0
		iindex = 0
		treatments = gene.getTreatments()
		# set union
		self.treatments = self.treatments | treatments
		if gene.strand == -1:
			# shallow copy of reversed list
			fbins = fbins[::-1]
			count = 0
			position = 0
			# reverse labelled bin positions, for validation
			for bin in fbins:
				bino = bin[1]
				bino.lower = position
				bino.upper = position + 1
				if count % 2 == 1:
					position += 1
				count += 1
		fbins, lastebins, lastibins = self.trimEIBins(fbins)
		for bin in fbins:
			(b_type, subbinlist) = bin
			try:
				if b_type == "intron":
					dropbinlist = self.intronbins[iindex]
					iindex += 1
				elif b_type == "exon" or b_type == "cexon":
					dropbinlist = self.exonbins[eindex]
					eindex += 1
				else:
					print "unknown bin type", b_type
					assert(False)
			except:
				settings.breakpoint()
			self.mergeSubBins(dropbinlist, subbinlist, gene.strand)
		self.mergeSubBins(self.lastexonbins, lastebins, gene.strand)
		self.mergeSubBins(self.lastintronbins, lastibins, gene.strand)

		
