#!/usr/bin/env python

import cintervaltree

pbins = [ (0,1), (10,11), (12,13), (14,15) ]

itree = cintervaltree.IntervalTree()
for pbin in pbins:
	itree.insert_interval(cintervaltree.Interval(pbin[0], pbin[1]))
	
print itree.find(8,15)
	

