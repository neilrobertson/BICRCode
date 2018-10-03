#!/usr/bin/env python

gene_ids = "gene_id_list.txt"
human_file = "human_all.txt"


genes_on = [ "AC116366.4", "IL4", "AC004775.1", "IL3", "AC034228.4", "AC063976.7", "SLC22A4", "SLC22A5", "AC063976.4", "RAD50", "PDLIM4", "UQCRQ", "IRF1", "AP006216.10", "ZNF259", "AP000936.1", "OR7E23P", "AP000279.69", "IL10RB", "IFNAR2", "IFNGR2", "IFNAR1", "AP000569.8", "C21orf119", "AP000297.1", "AP000271.1", "SON", "TMEM50B", "DONSON", "GART", "C21orf59", "RP11-143H17.1", "U52112.12", "IKBKG", "AVPR2", "CXorf12", "FUNDC2", "ATP6AP1", "BRCC3", "RPL10", "PLXNA3", "XX-FW81657B9.5", "TAZ", "EMD", "ATF4P", "FAM50A", "DKC1", "AF277315.13", "RENBP", "HCFC1", "FAM3A", "LAGE3", "SLC10A3", "MTCP1", "ARD1", "CTAG1B", "CTAG2", "FLNA", "MPP1", "RP11-115M6.4", "OR52D1", "TRIM22", "TRIM34", "TRIM6-TRIM34", "OR51B6", "OR51J1", "TRIM6", "AC104389.32", "AC104389.19", "HBB", "HBE1", "HBG1", "HBG2", "EVX1", "AC004080.15", "AC004079.7", "AC004079.4", "AC073472.1", "AC004079.2", "HOXA10", "TNNI2", "AC068580.6", "LSP1", "MRPL23", "AC139143.1", "CTSD", "AC051649.10", "TSPAN32", "RP11-545E17.13", "PPP2R4", "DOLPP1", "NUP188", "FAM73B", "TMEM15", "RP11-247A12.6", "SERF2", "HYPK", "PDIA3", "MAP1A", "AC023356.3", "MFAP1", "AP001187.11", "RASGRP2", "MEN1", "AP001187.10", "EHD1", "SF1", "RP3-477O4.13", "MMP24", "CEP250", "SPAG4", "RPL37P1", "ERGIC3", "C20orf52", "FER1L4", "GDF5", "NFS1", "RBM12", "CPNE1", "ITGB4BP", "FOXP4", "RP4-696P19.2", "RP11-298J23.7", "C6orf49", "FRS3" ]
genes_high = [ "RAD50", "PDLIM4", "UQCRQ", "AP000936.1", "IFNAR1", "AP000569.8", "C21orf119", "AP000297.1", "AP000271.1", "SON", "TMEM50B", "DONSON", "GART", "C21orf59", "BRCC3", "RPL10", "PLXNA3", "XX-FW81657B9.5", "TAZ", "EMD", "ATF4P", "FAM50A", "DKC1", "AF277315.13", "CTAG1B", "CTAG2", "FLNA", "MPP1", "RP11-115M6.4", "TRIM6", "AC104389.32", "HBB", "HBE1", "HBG1", "HBG2", "AC004079.7", "AC004079.4", "AC073472.1", "AC004079.2", "MRPL23", "AC139143.1", "PPP2R4", "DOLPP1", "NUP188", "FAM73B", "SERF2", "HYPK", "PDIA3", "MAP1A", "AC023356.3", "MFAP1", "AP001187.11", "AP001187.10", "EHD1", "SF1", "CEP250", "SPAG4", "RPL37P1", "ERGIC3", "C20orf52", "ITGB4BP", "FOXP4", "RP4-696P19.2", "RP11-298J23.7", "C6orf49" ]
genes_low = [ "AC116366.4", "IL4", "AC004775.1", "IL3", "AC034228.4", "AC063976.7", "SLC22A4", "SLC22A5", "AC063976.4", "IRF1", "AP006216.10", "ZNF259", "OR7E23P", "AP000279.69", "IL10RB", "IFNAR2", "IFNGR2", "RP11-143H17.1", "U52112.12", "IKBKG", "AVPR2", "CXorf12", "FUNDC2", "ATP6AP1", "RENBP", "HCFC1", "FAM3A", "LAGE3", "SLC10A3", "MTCP1", "ARD1", "OR52D1", "TRIM22", "TRIM34", "TRIM6-TRIM34", "OR51B6", "OR51J1", "AC104389.19", "EVX1", "AC004080.15", "HOXA10", "TNNI2", "AC068580.6", "LSP1", "CTSD", "AC051649.10", "TSPAN32", "RP11-545E17.13", "TMEM15", "RP11-247A12.6", "RASGRP2", "MEN1", "RP3-477O4.13", "MMP24", "FER1L4", "GDF5", "NFS1", "RBM12", "CPNE1", "FRS3" ]
genes_off = [ "IL5", "KIF3A", "AFF4", "APOC3", "C21orf55", "SYNJ1", "OLIG2", "C21orf62", "ITSN1", "AP000288.2", "OLIG1", "AP000269.3", "AP000569.2", "AP000282.3", "C21orf63", "H2AFB1", "CXorf2B", "AC092402.4", "CXorf2", "OPN1LW", "OPN1MW", "AC092402.5", "CTAG1A", "UBQLN3", "OR52A1", "OR52E1P", "OR51N1P", "OR51L1", "OR51A7", "OR52J3", "OR52J2P", "OR51F2", "OR51T1", "OR51H2P", "OR52U1P", "OR52J1P", "AC104389.16", "OR56B1", "OR52E3P", "HOXA9", "HOXA2", "HOXA7", "HOXA1", "HOXA6", "HOXA5", "TH", "IGF2", "INS", "IGF2AS", "AC051649.12", "C9orf106", "RP11-344B5.4", "RP11-344B5.3", "CKMT1A", "CKMT1B", "CKMT1", "AC011330.12", "NRXN2", "AP001092.5", "AP006288.1", "AP003774.5", "AP005273.1", "NCR2", "RP11-328M4.3", "RP11-298J23.5" ]

def processHumanData(filename):
	linesbyid = {}
	fh = open(filename, "r")
	header = fh.readline()
	for line in fh:
		line_sp = line.split("\t")
		this_id = line_sp[6]
		# this will replace lines, but any line will do for each id
		linesbyid[this_id] = line
#		this_list = linesbyid.setdefault(this_id, [])
#		this_list.append(line)
	fh.close()
	return linesbyid

def getIDs(filename):
	origids = {}
	fh = open(filename, "r")
	for line in fh:
		this_id = line[:-1]
		origids[this_id] = 1
	fh.close()
	return origids
	
ids = getIDs(gene_ids)
humandata = processHumanData(human_file)
missing = []
mfh = open("missing.txt", "w")
dfh = open("lookup_genes.txt", "w")
for this_id in ids:
	if this_id in humandata:
		line = humandata[this_id]
		if this_id in genes_high:
			activation = "high"
		elif this_id in genes_low:
			activation = "low"
		elif this_id in genes_on:
			activation = "on"
		elif this_id in genes_off:
			activation = "off"
		else:
			activation = "ambiguous"
#		for line in humandata[this_id]:
		print >> dfh, line[:-1] + "\t%s" % (activation)
	else:
		print >> mfh, this_id
		
dfh.close()
mfh.close()
