'''
Created on 17 Jun 2013

@author: mcbryan
'''

import getopt
import sys
from collections import defaultdict
from genemapping.chrmEnds import ChromosomeEnds
import csv

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["joinedfile=","binsize=","chra=","chrb=","genome=","minsize="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    chra = None
    chrb = None
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
        elif o=="--chra":
            chra = a
            print "Chr A: ", chra
        elif o=="--chrb":
            chrb = a
            print "Chr B: ", chrb
            #assert False, "Two different chromosomes not implemented yet"
        elif o=="--genome":
            genome = a
            print "Genome: ", genome
        elif o=="--minsize":
            minsize = int(a)
            print "Minsize: ", minsize
    
    assert chra != None
    assert joinedfile != None
    
    if chrb == None:
        chrb = chra
    
    chrmEnds = ChromosomeEnds(genome)
    
    matrix = defaultdict(int)
    
    reads = csv.reader(open(joinedfile, "r"), delimiter="\t")
    matrixcsv = csv.writer(open(joinedfile+"."+chra+".vs."+chrb+"."+str(binsize)+".matrix","w"),delimiter='\t')
    
    for read in reads:
        #print read[lread_chr],read[lread_pos],read[rread_chr],read[rread_pos]
                      
        lread_chr = read[lread_chr_col]
        rread_chr = read[rread_chr_col]
        lread_pos = int(read[lread_pos_col])
        rread_pos = int(read[rread_pos_col])
        
        lbin = lread_pos/binsize
        rbin = rread_pos/binsize
        
        # matrix is chrA vs chrB
        
        # on same chromosome, should be symmetric
        if (chra == chrb == lread_chr == rread_chr):
            
            if minsize != None and abs(rread_pos - lread_pos) < minsize:
                continue # skip anything less than a certain size
            
            matrix[(lbin, rbin)] += 1
            if (lbin != rbin): # but only if we aren't incrementing the same bin twice
                matrix[(rbin, lbin)] += 1
        
        elif (chra == lread_chr and chrb == rread_chr):
            matrix[(lbin, rbin)] += 1
            
        elif (chrb == lread_chr and chra == rread_chr):
            matrix[(rbin, lbin)] += 1
        
    for chrAPos in xrange(chrmEnds[chra]/binsize):
        row = []
        for chrBPos in xrange(chrmEnds[chrb]/binsize):
            row.append(matrix[(chrAPos,chrBPos)])
        matrixcsv.writerow(row)