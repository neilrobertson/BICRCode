#!/usr/bin/env python

import sys

crispr_filename = "/mnt/bam01/projects/Peter/CRISPR_EpigeneticLibrary/Cas9_exome_site_distribution_12122013/hg19_exome_12122013.SP.data.txt"
genes_filename = "/mnt/bam01/projects/Peter/CRISPR_EpigeneticLibrary/Chromatin_Genes.Gene_Symbol.txt"
ouput_filename = "/mnt/bam01/projects/Peter/CRISPR_EpigeneticLibrary/CasFinder.Results/Chromatin_Genes_27_01_2015_PeterAdams.CasFinder.data.txt"

genes_map = []
output_lines = []
header = None

with open(genes_filename, "r") as genes:
	for line in genes:
		gene = line.strip()
		genes_map.append(gene)


with open(crispr_filename, "r") as crispr_results:

	for i, line in enumerate (crispr_results):
		if i == 0:
			header = line.strip()
			print header
		else:
			line_parts = line.strip().split("\t")
			gene_id = line_parts[2] #.split(":")[0]
			if gene_id in genes_map:
				output_lines.append(line.strip())



with open(ouput_filename, "w") as output:

	output.write(header + "\n")
	for line in output_lines:
		output.write(line + "\n")