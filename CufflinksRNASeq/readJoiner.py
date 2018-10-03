'''
Created on 21 Mar 2011

@author: mcbryan
'''

import sys
import getopt
import subprocess
import csv
from sam.SamFormat import SAMFile, PairedReader
from StringIO import StringIO
from datastructures.genomeintervaltree import intervalsOverlap
import collections

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["bamfile="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    bamfile = None
    
    for o, a in opts:
        if o=="--bamfile":
            bamfile = a
    
    assert bamfile != None
        
    def getSAMreads(chrm, start, stop, bamfile):
        p = subprocess.Popen(['samtools', 'view', bamfile, chrm+":"+str(start)+"-"+str(stop)],
                                                  stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE)
        result, err = p.communicate()
        if p.returncode != 0:
            raise IOError(err)
        return result
    
    cols = {'chrm' : 0,
            'source' : 1,
            'type' : 2,
            'start': 3,
            'stop': 4,
            'details': 8}
    
    class NoRead():
        def getSegments(self):
            return []
    
    from Data import mergeInts
    
    from datastructures.cintervaltree import Interval
    
    # chr6:26124373-26139344 # spliced est in the middle
    # chr6:27106073-27114619
    for chrm,start,stop in [("chr6",27106073,27114619)]:
        
        for x,y in PairedReader(SAMFile(StringIO(getSAMreads(chrm,start,stop,bamfile)))).getReadPairs():
            
            # if both reads map in multiple places then we skip it
            if x.getNumberHits() > 1 and (y == None or y.getNumberHits() > 1):
                continue
            
            if y == None:
                y = NoRead()
            
            segments = set()
            segments.update(x.getSegments())
            segments.update(y.getSegments())
            
            
            starts = [ start for start,stop in segments ]
            stops = [ stop for start,stop in segments ]
            
            intervals = [Interval(start, stop) for start,stop in segments]

            maxcoord = max(max(starts),max(stops))
            mincoord = min(min(starts),min(stops))
            
            intervals = mergeInts(intervals)
           
            blockStarts = ",".join([ str(interval.start-mincoord) for interval in intervals ])
            blockSizes = ",".join([ str(interval.end-interval.start) for interval in intervals ])
           
            bedLine = [chrm, str(mincoord), str(maxcoord), ".", "0", ".", str(mincoord), str(maxcoord), "0,0,0", str(len(intervals)), blockSizes, blockStarts]
            print "\t".join(bedLine)