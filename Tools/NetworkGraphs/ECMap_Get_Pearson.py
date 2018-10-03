'''
Created on 9 Apr 2014

@author: johncole
'''

from scipy.stats.stats import pearsonr

import numpy as np
import sys
import getopt
import time
import multiprocessing as mp



if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["infile=","outfile=","minPCC=","threads=","header="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    inFile = None
    outFile = None
    minPCC = None
    threads = 2
    delimiter = "\t"
    header = True

    for o, a in opts:
        if o== "--infile":
            inFile = open(a).readlines()
        if o== "--outfile":
            outFile = open(a,'w')
        if o== "--minPCC":
            minPCC = float(a)  
        if o== "--threads":
            threads = int(a)  
        if o== "--header":
            header = bool(a)
            
    def calculate_pearsons(line,p1):
        lineParts = line.rstrip().split(delimiter)
        geneB = lineParts[0]
        lineParts.remove(geneB)
        p2 = np.array(lineParts).astype(np.float)
        r,p =  pearsonr(p1,p2)
        print geneB, r, p
        return geneB,r, p
        
    pool = mp.Pool(processes=threads)
    #We get the total number of genes:
    totalGenes = len(inFile)
    startTime = time.time()
    #We convert into pearson and print only where there is a significant correlation:
    counter = 0
    for line in range(totalGenes-1):
        if header == True and line == 0: pass
        else:
            counter +=1
        
            if counter%50 == 0:
                print "We have completed processing gene number %s" % str(counter)
            
            lineParts = inFile[line].rstrip().split(delimiter)
            geneA = lineParts[0]
            lineParts.remove(geneA)
            p1 = np.array(lineParts).astype(np.float) 
              
            results = [pool.apply(calculate_pearsons, args=(inFile[x],p1,)) for x in range(line+1,totalGenes-1)]
            for result in results:
                geneB, r, p = result[0], result[1], result[2]
                
                if abs(r) > minPCC: 
                    outFile.write(geneA + '\t' + geneB + '\t' + str(r) + '\t' + str(p) + '\n')
    
    outFile.flush()
    outFile.close()
    endTime = time.time() - startTime
    print "Time taken: " + endTime
