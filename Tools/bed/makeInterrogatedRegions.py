'''
Created on 21 Sep 2010

This was to make a list of interrogated regions for useq.

i.e. Chromosomes minus gaps

@author: mcbryan
'''

from genemapping.gaps import ChromosomeGaps
from genemapping.chrmEnds import ChromosomeEnds
from genemapping.chrList import ChrList
import collections

gaps = ChromosomeGaps("hg18")
ends = ChromosomeEnds("hg18")
chrs = ChrList("hg18")

chrCounter = collections.defaultdict(int)
regions = collections.defaultdict(list)

if __name__ == '__main__':
        
    for (chr,start,stop) in gaps.individualgaps:
        if chr in chrs:
            
            #print (chr,chrCounter[chr],start)
            regions[chr].append((chr,chrCounter[chr],start))
            chrCounter[chr] = stop
    
    for chr in chrs:
        #print (chr,chrCounter[chr],ends[chr])
        if chrCounter[chr] != ends[chr]:
            regions[chr].append((chr,chrCounter[chr],ends[chr]))
    
        for (chr,start,stop) in regions[chr]:
            print chr + "\t" + str(start) + "\t" + str(stop)