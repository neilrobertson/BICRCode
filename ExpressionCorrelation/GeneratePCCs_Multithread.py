'''
Created on 29 Jan 2015

@author: johncole
'''

#Inputs a matrix of expression (or any other) values, with the first column as genes names (or IDs). There is no header row.

from scipy.stats.stats import pearsonr
import sys
import getopt
import numpy as np
import multiprocessing as mp


def pearson(arrayX,arrayY,geneX,geneY):
    r,p =  pearsonr(arrayX,arrayY)
    return geneX,geneY,r,p


if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["infile=","outfile=","minPCC=","cores="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    inFile = None
    outFile = None
    minPCC = None
    cores = 1
    dataMatrix = []
    namesList = []

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
        if o== "--cores":
            cores = int(a)
            print "Number of threads = " + a
            
            
    pool = mp.Pool(processes=cores)
            
    print "#### Processing Data ####"
    
    #We get the total number of genes:
    totalGenes = len(inFile)
    print "Total genes in input:", totalGenes
    
    #We convert the infile into a numpy array:
    for line in inFile:
        tokens = line.rstrip().split('\t')
        dataMatrix.append(np.array(tokens[1:len(tokens)]).astype(np.float))
        namesList.append(tokens[0])
    print "Converted input file to data matrix"
    
    #We convert into pearson and print only where there is a correlation above the min PCC:
    counter = 0
        
    for xIndex in range(totalGenes-1):
        counter +=1 
    
        if counter%100 == 0:
            print "Genes Processed:", counter, "(" +str((float(counter)/float(totalGenes)*100)) + "%)"
        
        geneX = namesList[xIndex]
        arrayX = dataMatrix[xIndex]
        
        results = [pool.apply(pearson,args=(arrayX,dataMatrix[yIndex],geneX,namesList[yIndex],)) for yIndex in range(xIndex+1,totalGenes-1)]
        
        for result in results:
            if abs(result[2]) > minPCC: 
                outFile.write(str(result[0]) + '\t' + str(result[1]) + '\t' + str(result[2]) + '\t' + str(result[3]) + '\n')
       
    outFile.flush()
    outFile.close()
    
    print "#### Finished ####"
