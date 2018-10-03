#!/usr/bin/env python

import csv
from sets import Set

### Related and useful functionality in crosslookups.py

gene_coords = "/home/pzs/histone/ourmousegene.csv"
gene_activations = "/home/pzs/histone/act_list_new.txt"

reader = csv.reader(open(gene_coords, "rb"))

allgenes = Set()
for row in reader:
	allgenes.add(row[4])
	


reader = csv.reader(open(gene_activations, "rb"))

filedgenes = Set()
for row in reader:
	filedgenes.add(row[0])
	
diffrows = [ (gene, "unfiled") for gene in allgenes.difference(filedgenes) ]


writer = csv.writer(open("unfiled.txt", "wb"))
writer.writerows(diffrows)
	

