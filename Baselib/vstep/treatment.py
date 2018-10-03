from datastructures.cintervaltree import *
import csv

from datastructures.genomeintervaltree import GenomeIntervalTree

class Vstep(GenomeIntervalTree):
    def __init__(self, treatmentFileName, basewidth,limit=0):
        super(Vstep, self).__init__()
        self.basewidth = basewidth
        self.loadValues(treatmentFileName,limit)

    # get the values from the treatment
    def loadValues(self, valuesFile, limit = 0):
        values = {} # dictionary of chr -> interval tree (start, end, value)
        count = 0
        valuesInputFile = csv.reader(open(valuesFile, "r"), delimiter='\t')
        for row in valuesInputFile:
            #do int() casting when reading -> get it done early on
            count = count + 1
            chr = row[0]
            start = int(row[1])
            value = int(row[2])
            self.insertInterval(chr, start, start+self.basewidth, value)
            if limit > 0 and count > limit: # for debugging quickness
                break
        print "Values filed:", count
        return values
