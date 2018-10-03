'''
Created on 29 Jan 2015

@author: johncole
'''

import sys
import getopt

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["geneList=","outfile=","infile="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    inFile = None
    outFile = None
    geneLists = []

    for o, a in opts:
        if o== "--infile":
            inFile = open(a).readlines()
            print "File opened: " + a
        if o== "--outfile":
            outFile = open(a,'w')
            print "Outputting to: " + a
        if o== "--geneList":
            
            temp = {}
            tempInFile = open(a).readlines()
            
            for line in tempInFile:
                temp[line.rstrip()] = True
            geneLists.append(temp)
            print "Gene list added: " + a               

    
    print "#### Reading infile ####"
    
    for line in inFile:
        tokens = line.rstrip().split('\t')
        geneA = str(tokens[0])
        counter = 0
        for geneList in geneLists:
            if geneA in geneList.keys():
                outFile.write(geneA + '\t' + str(counter) + '\n')
            counter += 1

    print "Finished!"
