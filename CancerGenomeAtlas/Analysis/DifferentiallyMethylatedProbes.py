'''
Created on 28 Sep 2015

@author: neilrobertson
'''

import numpy as np
import scipy.stats as stats
from Filter import Filter
import sys, getopt

def getVariance(obs):
    p = 1.0
    try:
        chi2,p,dof,expected = stats.chi2_contingency(obs) #@UnusedVariable
    except ValueError:
        p = 1.0
    return p

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["infile1=","infile2=","outfile="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    filename_one = None
    filename_two = None
    outputFilename = None

    for o, a in opts:
        if o== "--infile1":
            filename_one = a
            print "File1: {0}".format(a)
        if o== "--infile2":
            filename_two = a
            print "File2: {0}".format(a)
        if o== "--outfile":
            outputFilename = a
            print "OutputFile: {0}".format(a)
            
assert filename_one
assert filename_two
assert outputFilename


DELIMITER = "\t"

cpgLoci = []

print "Buiding dictionary of first file data"
sample1Dict = {}
with open(filename_one, "r") as file_one:
    for i, line in enumerate(file_one):
        if i == 0: pass
        else:
            if i % 1000:
                print "Working on line {0}".format(str(i))
            lineParts = line.strip().split(DELIMITER)
            cpgLocus = lineParts[0]
            
            data = lineParts[1:] 
            data = Filter.convert_NAtoNaN(data)
            data = np.array(data).astype(np.float)   
            data = data[~np.isnan(data)]
            
            cpgLoci.append(cpgLocus)
            sample1Dict[cpgLocus] = data
            
print "Buiding dictionary of second file data"
sample2Dict = {}
with open(filename_two, "r") as file_two:
    for i, line in enumerate(file_two):
        if i == 0: pass
        else:
            if i % 1000:
                print "Working on line {0}".format(str(i))
            lineParts = line.strip().split(DELIMITER)
            cpgLocus = lineParts[0]
            
            data = lineParts[1:] 
            data = Filter.convert_NAtoNaN(data)
            data = np.array(data).astype(np.float)   
            data = data[~np.isnan(data)]
            
            sample2Dict[cpgLocus] = data
            
            
assert len(sample1Dict.keys()) == len(sample2Dict.keys())            
print "Built details from {0} cpg loci".format(str(len(sample1Dict.keys())))
            
with open(outputFilename, "w") as outputFile:
    
    headerLine = ["CpGRef", "median_one", "median_two", "geom_mean_1", "geom_mean_2", "median_diff", "mean_diff", "std_dev_1", "std_dev_2", "T", "mannu_PValue", "KS", "kolmog_PValue", "Welchs_T_PValue"]
    outputFile.write(DELIMITER.join(headerLine).strip() + "\n")
    
    for loci in cpgLoci:
        
        set_one = sample1Dict[loci]
        set_two = sample2Dict[loci]

        median_one = np.median(set_one)
        median_two = np.median(set_two)
        median_diff = median_two - median_one
        
        mean_one = stats.gmean(set_one)
        mean_two = stats.gmean(set_two)
        meanDiff = mean_two - mean_one
        
        std_one = np.std(set_one)
        std_two = np.std(set_two)
        
        KS, KOLMOGOROV_pValue = stats.ks_2samp(set_one, set_two)
        
        T, MANNWHITNEY_pValue = stats.mannwhitneyu(set_one, set_two)
        
        T_stat, ttest_PValue = stats.ttest_ind(set_one, set_two, equal_var=False)
        
        lociID = loci.strip().split("|")[0]     
        
        outputLine = [lociID, str(median_one), str(median_two), str(mean_one), str(mean_two), str(median_diff), str(meanDiff), str(std_one), str(std_two), str(T), str(MANNWHITNEY_pValue), str(KS), str(KOLMOGOROV_pValue), str(ttest_PValue)]
        outputFile.write(DELIMITER.join(outputLine).strip() + "\n")
        
print "complete!"