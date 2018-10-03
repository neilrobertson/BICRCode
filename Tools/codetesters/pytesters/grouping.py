#!/usr/bin/env python

from sets import Set

activations = {	"low"	: Set([1,2,3,4]),
				"high"	: Set([5,6]),
				"off"	: Set([7,8,9])}
				
grouplist = ["low", "high" ]
groupname = "on"

newact = Set()
for group in grouplist:
	acts = activations[group]
	newact = newact.union(acts)
activations[groupname] = newact

		
print activations
