'''
Created on 9 Jul 2013

@author: mcbryan
'''

import gzip
import mimetypes
import re

class FastAEntry(object):
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        assert header.startswith(">"), header + " : " + sequence
    
    def __str__(self):
        return self.header + "\n" + self.sequence
        
class FastAFile(object):

    def __init__(self, filename):
        self.filename = filename
        
    def readlines(self):
        if self.filename.endswith("gz") or "gzip" in mimetypes.guess_type(self.filename):
            self.reader = gzip.open(self.filename)
        else:
            self.reader = open(self.filename)
        
        header = None
        sequenceFragments = []
        
        for line in self.reader:
            
            if header == None:
                header = line.strip()
                continue
            
            if ">" in line:
                # new entry
                # output old one
                sequence = "".join(sequenceFragments)
            
                yield FastAEntry(header, sequence)
                
                # new entry uses the new header and empty the current list of sequence fragments
                header = line.strip()
                sequenceFragments = []
                
            else:
                sequenceFragments.append(line.strip())
            
        if len(sequenceFragments) > 0:
            yield FastAEntry(header, sequence)
        
        self.reader.close()
    
    def __iter__(self):
        return self.readlines()        
        
if __name__ == "__main__":
    
    for i in FastAFile("/mnt/50tb/privatedata/workspace/Baselib/fasta/TestFasta.fq"):
        print str(i)