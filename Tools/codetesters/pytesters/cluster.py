#!/usr/bin/env python

def findclusters(inlist, threshold):
	allclusters = []
	cluster = []
	adding = False
	i = 0
	while i < len(inlist)-1:
		thisval = inlist[i]
		nextval = inlist[i+1]
		if abs(nextval - thisval) < threshold:
			if not adding:
				cluster.append(thisval)
				adding = True
			cluster.append(nextval)
		else:
			if cluster:
				allclusters.append(cluster)
				cluster = []
			adding = False
		i += 1
	return allclusters
			

mylist = [1, 2, 4, 20, 25, 45, 47, 60, 66 ]
print findclusters(mylist, 6)
