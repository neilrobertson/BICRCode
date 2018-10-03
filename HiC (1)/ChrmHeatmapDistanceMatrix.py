'''
Created on 17 Jun 2013

@author: mcbryan
'''

import getopt
import sys
from collections import defaultdict
from genemapping.chrmEnds import ChromosomeEnds
from genemapping.chrList import ChrList
import csv

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["joinedfile=","genome=","minsize="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    joinedfile = None
    
    lread_chr_col = 3-1
    rread_chr_col = 8-1
    
    lread_pos_col = 4-1
    rread_pos_col = 9-1
    
    genome = "hg18"
    
    noSameChrm = True
    
    minsize = None
    
    for o, a in opts:
        if o=="--joinedfile":
            joinedfile = a
            print "Joined file: ", a
        elif o=="--genome":
            genome = a
            print "Genome: ", genome
        elif o=="--minsize":
            minsize = int(a)
            print "Minsize: ", minsize
    
    assert joinedfile != None
    
    chrmEnds = ChromosomeEnds(genome)
    chrmList = ChrList(genome)
    
    matrix = defaultdict(int)
    
    reads = csv.reader(open(joinedfile, "r"), delimiter="\t")
    matrixcsv = csv.writer(open(joinedfile+".interchromosomal.matrix","w"),delimiter='\t')
    
    for read in reads:
        #print read[lread_chr],read[lread_pos],read[rread_chr],read[rread_pos]              
        lread_chr = read[lread_chr_col]
        rread_chr = read[rread_chr_col]
        
        if noSameChrm and lread_chr == rread_chr:
            continue
        
        if minsize != None and lread_chr == rread_chr and abs(int(read[lread_pos_col]) - int(read[rread_pos_col])) < minsize:
            continue
               
        # should be symmetrical
        matrix[(lread_chr, rread_chr)] += 1
        if (lread_chr != rread_chr): # but only if we aren't incrementing the same bin twice
            matrix[(rread_chr, lread_chr)] += 1
    
    def geomean(nums):
        return (reduce(lambda x, y: float(x)*float(y), nums))**(1.0/len(nums))
    
    chrmList.remove("chrM")
        
    for chrA in chrmList:
        row = []
        for chrB in chrmList:
            row.append(float(matrix[(chrA,chrB)]))
                       #/(chrmEnds[chrA]*chrmEnds[chrB]))
        matrixcsv.writerow(row)
        print "Goood grief, im done!"