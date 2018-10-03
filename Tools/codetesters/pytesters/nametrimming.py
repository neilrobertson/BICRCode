#!/usr/bin/env python

mystr = "2344_36_K79me1_0--2344_38_K79me1_5--2344_37_K79m"
mystr2 = "2344_3_RabbitIgG_10--2344_28_RabbitIgG_10--2344_4_RabbitI"


def trimFilename(filename, length):
	trimstrings = ["-", "_"]
	trimtos = []
	for trimstring in trimstrings:
		thistrimto = filename.find(trimstring, length)
		trimtos.append(thistrimto)
	trimto = min(trimtos)
	return filename[:trimto]

print len(mystr)
print len(mystr2)
print trimFilename(mystr, 30)

