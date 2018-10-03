#!/usr/bin/env python

# num is the number of bins. Half will be exon and half intron
def makeExonIntronXAxisLabels(num):
	reverse = False
	labels = []
	# we only go from 1 to num because we're going to add two more for the last entries
	# that is, if we're not doing a reverse plot...
	if reverse:
		num += 1
	pairs = (num / 2) + 1
	for i in range(1, pairs):
		if reverse:
			elabel = "e -" + str(i)
			ilabel = "i -" + str(i)
		else:
			elabel = "e" + str(i)
			ilabel = "i" + str(i)
		labels.extend([elabel, ilabel])
	if num % 2 == 1:
		labels = labels[:-1]
	if not reverse:
		labels.extend([ "Ilast", "Elast"])
#	assert (len(labels) == num + 2)
	return labels

def makeExonIntronXAxisLabels2(num):
	labels = []
	reverse = True
	for i in range(num):
		number = str((i/2)+1)
		if reverse:
			number = " -" + number
		if i % 2 == 0:
			label = "e" + number
		else:
			label = "i" + number
		labels.append(label)
	if not reverse:
		labels.extend([ "Ilast", "Elast"])
	return labels	

print makeExonIntronXAxisLabels2(9)
