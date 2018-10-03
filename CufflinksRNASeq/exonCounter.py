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
        opts, args = getopt.getopt(sys.argv[1:], "", ["exonsdesc=","bamfile="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    exonsdesc = None
    bamfile = None
    
    for o, a in opts:
        if o=="--exonsdesc":
            exonsdesc = a
        elif o=="--bamfile":
            bamfile = a
    
    assert exonsdesc != None
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
    
    for row in csv.reader(open(exonsdesc,"r"),delimiter="\t"):
        
        chrm, type, start, stop, details = row[cols['chrm']], row[cols['type']], int(row[cols['start']]), int(row[cols['stop']]), row[cols['details']]
        
        if type == "exonic_part":
            
            fragments = 0
            
            for x,y in PairedReader(SAMFile(StringIO(getSAMreads(chrm,start,stop,bamfile)))).getReadPairs():
                
                # if both reads map in multiple places then we skip it
                if x.getNumberHits() > 1 and (y == None or y.getNumberHits() > 1):
                    continue
                
                if y == None:
                    y = NoRead()
                
                segments = set()
                segments.update(x.getSegments())
                segments.update(y.getSegments())
                
                for segmentstart,segmentend in segments:
                    if intervalsOverlap(start,stop,segmentstart,segmentend):
                        fragments += 1
                        break
            
            records = details.split(";")
            attributes = {}
            for record in records:
                key, value = record.strip().split(" ")
                attributes[key] = value.replace("\"","")
            
            print attributes['gene_id'], attributes['exonic_part_number'], stop-start, fragments
    