'''
Created on 20 Mar 2014

@author: johncole
'''
import scipy.stats

infile = open("/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/CpG_context_Pooled.binom.0.01FDR_CoverageCorrected.txt").readlines()
outfile = open("/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/CpG_context_Pooled.binom.0.01FDR_CoverageCorrected_CpGDiff.csv","w")

for line in infile:
    tokens = line .rstrip().split('\t')
    
    if int(tokens[2]) > 9 and int(tokens[3]) > 9 and int(tokens[4]) > 9 and int(tokens[5]) > 9:
        oddsratio, pvalue = scipy.stats.fisher_exact([[int(tokens[2]), int(tokens[3])], [int(tokens[4]), int(tokens[5])]])
        outfile.write(line.rstrip() + "\t" + str(oddsratio) + "\t" + str(pvalue) + "\n")
        
outfile.flush()
outfile.close()