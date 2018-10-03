'''
Created on 29 Jan 2015

@author: johncole
'''

#Inputs a matrix of expression (or any other) values, with the first column as genes names (or IDs). There is no header row.

from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
import sys
import getopt
import numpy as np

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["infile=","outfile=","minPCC=", "method="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    inFile = None
    outFile = None
    minPCC = None
    method = "both"
    dataMatrix = {}

    for o, a in opts:
        if o== "--infile":
            inFile = open(a).readlines()
            print "File opened: " + a
        if o== "--outfile":
            outFile = open(a,'w')
            print "Outputting to: " + a
        if o== "--minPCC":
            minPCC = float(a)    
            print "Printing all correlations >=" + a
        if o== "--method":
            method = a
            assert a.lower() in ["p","s","b","pearson", "spearman", "both"]
            print "Printing results based on:" + a   
            
    print "#### Processing Input Data ####"
    
    #We get the total number of genes:
    totalGenes = len(inFile)
    print "Total genes in input:", totalGenes
    print "Converting input"
    
    #We convert the data into numpy arrays:
    for line in inFile:
        tokens = line.rstrip().split('\t')
        gene = tokens[0]
        tokens.remove(gene)
        
        try:
            dataMatrix[gene] = np.array(tokens).astype(np.float)
        except:
            print "Data error: ",gene, len(tokens), tokens

    print "#### Calculating PCCs ####"
  
    #We calculate the PCCs and print only where above the threshold
    originalKeys = dataMatrix.keys()
    
    #We get gene A:
    for keyX in originalKeys:
        if keyX.strip() == "DNMT3B.cnv":
            print "YA DANCER!!"
            geneX = dataMatrix[keyX]
            del dataMatrix[keyX]
            
            #We provede a count of the remaining genes:
            if len(dataMatrix.keys()) % 100 == 0:
                print "Genes remaining = ",len(dataMatrix.keys())
            
            #We get keyY
            for keyY in dataMatrix.keys():
                geneY = dataMatrix[keyY]
                
                if method == "p" or method == "pearson":
                    r1,p1 =  pearsonr(geneX,geneY)
                    if abs(r1) > minPCC: 
                        outFile.write(keyX + '\t' + keyY + '\t' + str(r1) + '\t' + str(p1) + '\n')
                 
                elif method == "s" or method == "spearman":
                    r2,p2 =  spearmanr(geneX,geneY)
                    if abs(r2) > minPCC: 
                        outFile.write(keyX + '\t' + keyY + '\t' + str(r2) + '\t' + str(p2) + '\n')               
                else:
                    r1,p1 =  pearsonr(geneX,geneY)
                    r2,p2 =  spearmanr(geneX,geneY)            
                    if abs(r1) > minPCC and abs(r2) > minPCC: 
                        outFile.write(keyX + '\t' + keyY + '\t' + str(r2) + '\t' + str(p2) + '\n')            
            
    outFile.flush()
    outFile.close()
    
    print "#### Finished ####"
