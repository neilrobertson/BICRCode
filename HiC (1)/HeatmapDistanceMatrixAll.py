'''
Created on 17 Jun 2013

@author: mcbryan
'''

import getopt
import sys
from collections import defaultdict
from genemapping.chrmEnds import ChromosomeEnds
import csv
from genemapping.chrList import ChrList

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["joinedfile=","binsize=","genome=","minsize="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    binsize = 1000000
    joinedfile = None
    
    lread_chr_col = 3-1
    rread_chr_col = 8-1
    
    lread_pos_col = 4-1
    rread_pos_col = 9-1
    
    genome = "hg18"
    
    minsize = None
    
    for o, a in opts:
        if o=="--joinedfile":
            joinedfile = a
            print "Joined file: ", a
        elif o=="--binsize":
            binsize = int(a)
            print "Binsize: ", a
            #assert False, "Two different chromosomes not implemented yet"
        elif o=="--genome":
            genome = a
            print "Genome: ", genome
        elif o=="--minsize":
            minsize = int(a)
            print "Minsize: ", minsize
    
    assert joinedfile != None
    
    chrmEnds = ChromosomeEnds(genome)
    chrmList = ChrList(genome)
    
    chrmMatrix = defaultdict(dict)
    for chrmA in chrmList:
        for chrmB in chrmList:
            chrmMatrix[chrmA][chrmB] = defaultdict(int)
            chrmMatrix[chrmB][chrmA] = defaultdict(int)
    
    reads = csv.reader(open(joinedfile, "r"), delimiter="\t")
    
    for read in reads:
        #print read[lread_chr],read[lread_pos],read[rread_chr],read[rread_pos]
        
        lread_chr = read[lread_chr_col]
        rread_chr = read[rread_chr_col]
        lread_pos = int(read[lread_pos_col])
        rread_pos = int(read[rread_pos_col])
        
        if lread_chr == "*" or rread_chr == "*":
            continue
        
        lbin = lread_pos/binsize
        rbin = rread_pos/binsize
        
        # matrix is chrA vs chrB
        matrix = chrmMatrix[lread_chr][rread_chr]
        
        # on same chromosome, should be symmetric
        if (lread_chr == rread_chr):
            
            if minsize != None and abs(rread_pos - lread_pos) < minsize:
                continue # skip anything less than a certain size
            
            matrix[(lbin, rbin)] += 1
            if (lbin != rbin): # but only if we aren't incrementing the same bin twice
                matrix[(rbin, lbin)] += 1
                
            continue
        
        matrix[(lbin, rbin)] += 1
        
        # do the symmetric version too
        matrix = chrmMatrix[rread_chr][lread_chr]
        matrix[(rbin, lbin)] += 1
    
    chrmList.remove("chrM")
    
    for chra in chrmList:
        for chrb in chrmList:
            matrixcsv = csv.writer(open(joinedfile+"."+chra+".vs."+chrb+"."+str(binsize)+".matrix","w"),delimiter='\t')
            
            matrix = chrmMatrix[chra][chrb]
                
            for chrAPos in xrange(chrmEnds[chra]/binsize):
                row = []
                for chrBPos in xrange(chrmEnds[chrb]/binsize):
                    row.append(matrix[(chrAPos,chrBPos)])
                matrixcsv.writerow(row)
                
            