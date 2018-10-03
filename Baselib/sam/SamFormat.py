import csv

# A - char
# i - signed int
# f - float
# Z - string
# H - hex
def processOptionalField(field):
    (name, type, value) = field.split(":")
    if type == "i":
        value = int(value)
    elif type == "f":
        value = float(value)
    elif type == "H":
        value = int(value,16)
    return (name,value)

class SAMEntry():
    def __init__(self, entry):
        keys = {"key" : 0,
                "flags" : 1,
                "chrm" : 2,
                "start": 3,
                "mapping-quality-score":4,
                "cigar":5,
                "matechrm":6, # '*' if single ended
                "matepos":7, # 0 if singe ended
                "inferredinsert":8, # distance between pairs
                "sequence":9,
                "quality":10,
                # ... <anything additional to this is optional field identified by name>
                "optional-starts-at":11 
                }
        
        if entry[0].startswith("@"):
            # header entry
            self.header = entry
            return
        else:
            self.header = None
        
        self.key = entry[keys["key"]]
        self.flags = int(entry[keys["flags"]])
        self.chrm = entry[keys["chrm"]]
        self.start = int(entry[keys["start"]])-1 # SAM format is 1 based, convert to 0 based
        self.cigar = entry[keys["cigar"]]
        self.sequence = entry[keys["sequence"]]
        
        self.optionalFields = {}
        
        for i in range(keys["optional-starts-at"],len(entry)):
            name, value = processOptionalField(entry[i])
            self.optionalFields[name] = value
    
    def __repr__(self):
        if self.isHeader():
            return str(self.header)
        else:
            return self.key +" / " + str(self.getNumberHits()) + " / "+self.chrm + ":"+str(self.start) + " / " + self.cigar + " / " + str(self.getSegments())
    
    def checkFlag(self,bitmask):
        return self.flags & bitmask
    
    def isHeader(self):
        return self.header != None
    
    def getNumberHits(self):
        assert not self.isHeader()
        if "NH" in self.optionalFields:
            return self.optionalFields["NH"]
        else:
            return 1 # if marker absent assume one hit
    
    def getStrand(self):
        assert not self.isHeader()
        if "XS" in self.optionalFields:
            strand = self.optionalFields["XS"]
            assert strand == "+" or strand == "-" or strand == "."
            return strand
        else:
            return "."
    
    def getSegments(self):
        assert not self.isHeader()
        segments = []
        pos = self.start
        current = ''
        for char in self.cigar:
            # N = skip this number of bases, I = there's bases in the reference not in the read (so still skip this space)
            if char == 'N' or char == 'D':
                assert len(current) > 0, "N/D without number of bases?"
                pos += int(current)
                current = ''
            # insertion into the read (i.e. the read has extra bases that arent in the reference
            elif char == 'I':
                # dont update anything
                assert len(current) > 0, "I without number of bases?"
                continue
            # matching bases
            elif char == 'M':
                assert len(current) > 0, "M without number of bases?"
                length = int(current)
                segments.append((pos, pos + length))
                current = ''
                pos += length
            elif char in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']:
                # numerics
                current += char
            else:
                assert False, "I'm not sure how to deal with this CIGAR string" + self.cigar
        return segments 



class SAMFile():
    # assumes tophat style same files
    def __init__(self, f):
        if isinstance(f,basestring):
            self.reader = csv.reader(open(f), delimiter="\t")
        else:
            self.reader = csv.reader(f,delimiter="\t")
        
    def next(self):
        return SAMEntry(self.reader.next())

    def __iter__(self):
        return self


class PairedReader():
    def __init__(self,samfile):
        self.samfile = samfile
        self.unpaired = {}
    
    def getReadPairs(self):
        for nextread in self.samfile:
            if nextread.key in self.unpaired:
                one, two = nextread, self.unpaired[nextread.key]
                del self.unpaired[nextread.key]
                yield one,two
            else:
                self.unpaired[nextread.key] = nextread

        # no more reads to get from the sam file so go to unpaired reads
        for key in self.unpaired:
            yield self.unpaired[key], None
        
        # reset list of unpaired reads
        self.unpaired = {}
        
        # no more unpaired reads
        raise StopIteration
    
    
if __name__ == "__main__":

    for x in PairedReader(SAMFile("/mnt/50tb/privatedata/Taranjit/Ribodepleted-RNA-Seq/exon-counting/PD32.lane1.sam")).getReadPairs():
        print x
        #print one,two