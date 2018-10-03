'''
Created on 23 Aug 2010

@author: mcbryan
'''
import sys
import getopt
from csvfile.indexedcsv import IndexedCSV
import collections
from fisher import cfisher

def fisherExact(uphasTF,upnoTF,downhasTF,downnoTF):    
    p = cfisher.pvalue(uphasTF,upnoTF,downhasTF,downnoTF).two_tail
    #return -10.0*math.log10(p)
    return p

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["gene-expr-transcriptionfactor-file="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    for o, a in opts:
        if o=="--gene-expr-transcriptionfactor-file":
            infile = a
            print "Matrix Dir:",infile            
    assert infile!=None
    
    ensemblidcol = "ensemblid"
    
    inCSV = IndexedCSV(infile,keyPos=1)

    keys = inCSV.keys
    
    transcriptionFactors = keys[keys.index("significant")+1:]
    
    print transcriptionFactors
    
    assert len(set(transcriptionFactors))==len(transcriptionFactors)
    
    transcriptionFactorUp = collections.defaultdict(int)
    transcriptionFactorDown = collections.defaultdict(int)
    transcriptionFactorNC = collections.defaultdict(int)
    
    totalUp = 0
    totalDown = 0
    totalNC = 0
    
    # count how often the TFs appear
    for row in inCSV:
        
        significant = False
        if inCSV[row]['significant'] == 'yes':
            significant = True
        direction = "up" if float(inCSV[row]['log2(fold_change)'])>=0 else "down"
        
        if significant and direction=="up":
            totalUp += 1
        elif significant and direction=="down":
            totalDown += 1
        else:
            totalNC += 1
        
        for transcriptionFactor in transcriptionFactors:
            
            if significant and direction=="up" and int(inCSV[row][transcriptionFactor]) >=10:
                transcriptionFactorUp[transcriptionFactor] += min(int(inCSV[row][transcriptionFactor]),1)
            elif significant and direction=="down" and int(inCSV[row][transcriptionFactor]) >=10:
                transcriptionFactorDown[transcriptionFactor] += min(int(inCSV[row][transcriptionFactor]),1)
            else:
                transcriptionFactorNC[transcriptionFactor] += min(int(inCSV[row][transcriptionFactor]),1)
   
    
    for transcriptionFactor in transcriptionFactors:
        print transcriptionFactor, \
            transcriptionFactorNC[transcriptionFactor], \
            totalNC - transcriptionFactorNC[transcriptionFactor], \
            transcriptionFactorUp[transcriptionFactor], \
            totalUp - transcriptionFactorUp[transcriptionFactor], \
            transcriptionFactorDown[transcriptionFactor], \
            totalDown - transcriptionFactorDown[transcriptionFactor], \
            fisherExact(transcriptionFactorUp[transcriptionFactor],totalUp-transcriptionFactorUp[transcriptionFactor],transcriptionFactorDown[transcriptionFactor],totalDown-transcriptionFactorDown[transcriptionFactor]), \
            fisherExact(transcriptionFactorNC[transcriptionFactor], totalNC - transcriptionFactorNC[transcriptionFactor], transcriptionFactorUp[transcriptionFactor]+transcriptionFactorDown[transcriptionFactor], totalDown-transcriptionFactorDown[transcriptionFactor]+totalUp-transcriptionFactorUp[transcriptionFactor])
            