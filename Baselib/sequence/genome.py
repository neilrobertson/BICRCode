import csv
import os
import sys
import re
from genemapping.chrmEnds import ChromosomeEnds

# realises that there is a line break every 50 bases which uses up one byte in the input
def fixedPosition(position, build = "hg19"):
    numb_chars_per_line = 50
    if build.strip() == "mm9":
        numb_chars_per_line = 70
    return (((position/numb_chars_per_line)*(numb_chars_per_line+1)) + position%numb_chars_per_line)
    
class UnknownChromosomeException(Exception):
    def __init__(self,  chr):
        self.chr = chr
    def __str__(self):
        return self.chr + " not a known chromosome"


class Genome():
    def __init__(self, genomeBuild = "hg19"):        
        self.sequenceBase = os.path.expanduser("~/mount/publicdata/" + genomeBuild + "/assembly/fasta/")
        self.build = genomeBuild
    
    def getSequence(self, chr, *arguments):
        # (chr, (start,end))
        # or
        # (chr,[start,end])
        
        if len(arguments)==1:
            (start, end) = arguments[0]
        elif len(arguments)==2:
            start = arguments[0]
            end = arguments[1]
        else:
            assert False
            
        if start <0:
            start =0

        fastafile = self.sequenceBase + chr + ".fa"
        
        try:
            with open(fastafile, "rb") as chromosome:
                
                assert(start>=0), chr +":" + str(start) + "-" + str(end)
                assert(end-start>=0),"Start:"+str(start)+", End:"+str(end)
                
                chromosome.readline() # read in the header line which looks like ">chrBlah"
                newstart = fixedPosition(start, self.build)
                newend = fixedPosition(end, self.build)
                chromosome.seek(newstart, os.SEEK_CUR)            
                sequence = chromosome.read(newend-newstart).replace("\n", "")
                assert len(sequence) <= end-start
                return sequence
        except IOError:
            raise UnknownChromosomeException(chr)

    @staticmethod
    def IUPACtoRegex(iupacString):
        map = {        'A' : 'Aa',  # adenine
                       'C' : 'Cc',  # cytosine
                       'G' : 'Gg',  # guanine
                       'T' : 'Tt',  # thymine
                       'R' : 'GgAa',  # G A (purine)
                       'Y' : 'TtCc',  # T C (pyrimidine)
                       'K' : 'GgTt',  # G T (keto)
                       'M' : 'AaCc',  # A C (amino)
                       'S' : 'GgCc',  # G C (strong bonds)
                       'W' : 'AaTt',  # A T (weak bonds)
                       'B' : 'GgTtCc',  # G T C (all but A)
                       'D' : 'GgAaTt',  # G A T (all but C)
                       'H' : 'AaCcTt',  # A C T (all but G)
                       'V' : 'GgCcAa',  # G C A (all but T)
                       'N' : 'AaGgCcTt' # A G C T (any)
                       }
        # if the iupacString has anything that's not in the mapping it'll throw a keyerror
        regexStringArray = ['['+map[x]+']' for x in iupacString.upper()]
        return ''.join(regexStringArray)

    @staticmethod
    def reverseBrackets(character):
        mapping = {'[':']',
               ']':'['}
        return mapping[character] if character in map else character
        
    # this works for sequences and regex'es (hence the [] bit)
    @staticmethod
    def complement(string):
        mapping = {   'A':'T',
                      'T':'A', 
                      'G':'C', 
                      'C':'G', 
                      
                      'a':'t', 
                      't':'a', 
                      'g':'c', 
                      'c':'g',
                      
                      'R':'Y', #GA -> CT
                      'Y':'R', #CT -> GA
                      'K':'M', #GT -> CA
                      'M':'K', #CA -> GT
                      'S':'S', #GC -> CG
                      'W':'W', #AT -> TA                                 
                      
                      'r':'y', #GA -> CT
                      'y':'r', #CT -> GA
                      'k':'m', #GT -> CA
                      'm':'k', #CA -> GT
                      's':'s', #GC -> CG
                      'w':'w', #AT -> TA  
 
                      'N':'N',
                      'n':'n',
                      
                      '[':'[', 
                      ']':']'}
        # if the regexString has anything that's not in the mapping it'll throw a keyerror
        stringArray = [mapping[x] for x in string]
        return ''.join(stringArray)
    
    # reverse the string, flip any regex braces leaving all the bases unaltered
    @staticmethod
    def reverse(string):
        mapping = {   '[':']', 
                      ']':'['}
        stringArray = [mapping[x] if x in mapping else x for x in string[::-1]]
        return ''.join(stringArray)
    
    # assumes sequence is always ordered such that for coord of the sequence in the genome : coord(left) <= coord(right)
    # that is to say 5' -> 3' for the + strand
    @staticmethod
    def sequenceHasIUPAC(sequence, iupacString, strand=None):
        assert strand=='+' or strand=='-' or strand == None
        
        if strand == None:
            return (Genome.sequenceHasIUPAC(sequence,iupacString,strand='+')
                    or Genome.sequenceHasIUPAC(sequence,iupacString,strand='-'))
        
        if strand=='-':
            # reverse it and complement it
            sequence = Genome.complement(Genome.reverse(sequence))
        
        regex =  Genome.IUPACtoRegex(iupacString)
        if re.search(regex,  sequence) == None:
            return False
        else:
            return True


class MemoryGenome(Genome):
    def __init__(self, genomeBuild = "hg18"):
        self.sequenceBase = os.path.expanduser("~/mount/publicdata/" + genomeBuild + "/assembly/fasta/")
        
        self.chromosomeEnds = ChromosomeEnds(genomeBuild)
        
        self.chrms = {}
    
    def loadChrm(self,chr):
        fastafile = self.sequenceBase + chr + ".fa"
        try:
            with open(fastafile, "rb") as chromosome:
                
                chromosome.readline() # read in the header line which looks like ">chrBlah"
                sequence = chromosome.read().replace("\n", "")
                assert len(sequence)==self.chromosomeEnds[chr]
                self.chrms[chr]=sequence
        except IOError:
            raise UnknownChromosomeException(chr)
    
    def getSequence(self,chr, *arguments):
        # (chr, (start,end))
        # or
        # (chr,[start,end])
        
        if len(arguments)==1:
            (start, end) = arguments[0]
        elif len(arguments)==2:
            start = arguments[0]
            end = arguments[1]
        else:
            assert False
            
        if start <0:
            start =0
        
        assert(start>=0), chr +":" + str(start) + "-" + str(end)
        assert(end-start>=0),"Start:"+str(start)+", End:"+str(end)
        
        if chr not in self.chrms:
            self.loadChrm(chr)
        
        return self.chrms[chr][start:end]


if __name__ == "__main__":
    
    genome = Genome("hg18")
    memgenome = MemoryGenome("hg18")
    
    assert genome.getSequence("chrM", 0, 10) == memgenome.getSequence("chrM", 0, 10)
    assert genome.getSequence("chrM", 0, 20) == memgenome.getSequence("chrM", 0, 20)
    assert genome.getSequence("chrM", 40, 60) == memgenome.getSequence("chrM", 40, 60)
    assert genome.getSequence("chrM", 16070, 17650) == memgenome.getSequence("chrM", 16070, 17650)
    assert genome.getSequence("chr1", 16070, 17650) == memgenome.getSequence("chr1", 16070, 17650)
    assert genome.getSequence("chr1", 247249710, 247249730) == memgenome.getSequence("chr1", 247249710, 247249730)
    
    import timeit
        
    print timeit.timeit(stmt='genome.getSequence("chr1", 247249710, 247249730)', setup='from __main__ import genome', number=1000)
    
    print timeit.timeit(stmt='memgenome.getSequence("chr1", 247249710, 247249730)', setup='from __main__ import memgenome', number=1000)

    
    
    
