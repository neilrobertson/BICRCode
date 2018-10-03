#!/usr/bin/env python

import glob
import os

filenames = glob.glob("/home/pzs/histone/HISTONE_DATA/data/final/*.final.zscored.txt")
filenames = [ os.path.basename(f) for f in filenames ]


prefixes = [ "K562_", "U937_", "Monocyte_" ]

# takes a list of pointset "prefixes"
# combines pointsets that have the same name after removal of the prefix
# for example given a prefix list of ["K562_", "U937_" ], we should combine K562_H3K27Me1 and U937_H3K27Me1
def combineByPrefix(filenames, prefixes):
	# this will eventually have a base name (without any prefix)
	# hashed to a list of pointsetnames that have this basename with a prefix
	combinations = {}
	for psname in filenames:
		for prefix in prefixes:
			if psname.startswith(prefix):
				strippedname = psname[len(prefix):]
				comblist = combinations.setdefault(strippedname, [])
				comblist.append(psname)
	return combinations
	
combs = combineByPrefix(filenames, prefixes)

for combname in combs:
	print combname
	print combs[combname]
	print

