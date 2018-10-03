'''
Created on 20 Mar 2014

@author: johncole
'''
import sys, getopt
import scipy.stats

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:],"",["inputFile=", "outputFile=", "minCoverage="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    inFile = None
    outFile= None
    minCoverage = 10
    
    for o, a in opts:
        if o=="--inputFile":
            inFile = a
        elif o=="--outputFile":
            outFile = a
        elif o=="--minCoverage":
            minCoverage = int(a)
    
    assert inFile != None
    assert outFile != None

with open(inFile, "r") as inputFile:
    outputFile = open(outFile, "w")
    for line in inputFile:
        (chrm, position, lmeth, lunmeth, rmeth, runmeth) = line .rstrip().split('\t')
        
        if (int(lmeth) + int(lunmeth)) >= minCoverage and (int(rmeth) + int(runmeth)) >= minCoverage:
            oddsratio, pvalue = scipy.stats.fisher_exact([[int(lmeth), int(lunmeth)], [int(rmeth), int(runmeth)]])
            outputFile.write(line.rstrip() + "\t" + str(oddsratio) + "\t" + str(pvalue) + "\n")
            
    outputFile.flush()
    outputFile.close()