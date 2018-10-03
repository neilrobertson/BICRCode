'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv
from sam.SamFormat import SAMFile
from sequence.genome import MemoryGenome
from genemapping.chrmEnds import ChromosomeEnds

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["sam="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    debug = False
    
    infile = None
    strand = "+"
    
    for o, a in opts:
        if o=="--sam":
            infile = a
            print "SAM:", a
    
    assert infile != None
    
    with open(infile+".methcalls","w") as outfile:
        
        csvout = csv.writer(outfile,delimiter="\t")
    
        genome = MemoryGenome("hg18")
        ends = ChromosomeEnds("hg18")
        sam = SAMFile(infile)
        
        def findGenerator(haystack,needle):
            index = 0
            while True:
                index = haystack.find(needle,index)
                if index == -1:
                    break
                else:
                    yield index
                    index += 1
        
        for read in sam:
    
            # chromosomes we don't know about (chrL is probably lambda)
            if read.chrm not in ends:
                continue
            
            # check strand of the read
            strand = "-" if read.checkFlag(0x10) else "+"
    
            # get the genomic sequence (extended by 1 base on each side)
            # check the genomic sequence for CG dinucleotides
            # for each CG, extract the base from the sequence read corresponding to the C (make sure we do this right for both strands)
            #  if it is a C then it is methylated
            #  if it is a T then it is unmethylated
            # only valid for Lister et al. style sequencing (where complementary sequences to original DNA are not sequenced)
            
            # we want to get the genomic sequence + 1 base on each side 
            
            
            # assume we can't go off the beginning and the end of the sequence at the same time
            # ie. our chromosomes are comparatively large compared to our reads
            assert read.start-1 > 0 or read.start+len(read.sequence)+1 <= ends[read.chrm]
            
            if read.start-1 < 0:
                # we will go off the beginning of the chromosome, pad with one N
                genomeseq = "N" + genome.getSequence(read.chrm,read.start,read.start+len(read.sequence)+1)
                
            elif read.start+len(read.sequence)+1 > ends[read.chrm]:
                # we will go off the end of the chromosome, pad with one N
                genomeseq = genome.getSequence(read.chrm,read.start-1,read.start+len(read.sequence)) + "N"
            
            else:
                genomeseq = genome.getSequence(read.chrm,read.start-1,read.start+len(read.sequence)+1)
    
            # make the two sequences comparable in terms of character set (all uppercase) + start positions
            read.sequence = "N"+read.sequence.upper()+"N"
            genomeseq = genomeseq.upper()
    
            # do a check to see if there are any CG's in there first (slower than not checking of course)
            # only searches genomic forward strand but this is fine since CG's are the same on both strands
            if "CG" in genomeseq:
                
                if debug:
                    print
                    print read.chrm, read.start, len(read.sequence)+read.start
                    print read.sequence
                    print genomeseq
                    print "CG".join(["-"*len(seq) for seq in genomeseq.split("CG")])
                
                # outputs (C,G) locations
                locs = [(C,C+1) for C in findGenerator(genomeseq,"CG")]
                
                if strand == "+":
                    bases = [(read.start+C,read.sequence[C:C+1]) for C,G in locs]
                else:
                    # we want the G from the CG (which is a C on the opposite strand and which is one base along)
                    # note that the sequence in the SAM file is complemented compared to the genome. i.e. it's the actual
                    # sequence from the sequencer and will still have C or T as the basis for meth / unmeth calls
                    bases = [(read.start+G,read.sequence[G:G+1]) for C,G in locs]
                    
                for pos, base in bases:
                    if base in ["C","T"]: # ignore anything that's got an N or a SNP at that position
                        # we can make a meth call
                        # C methylated, T unmethylated
                        methCall = "z" if base == "T" else "Z"
                        methState = "-" if base == "T" else "+"
                        
                        csvout.writerow([read.key,methState,read.chrm,pos,methCall])