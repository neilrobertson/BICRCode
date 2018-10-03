import getopt
import sys
import csv
import math
from datetime import datetime


def parseMethLine(line):
    (chrm,coord,meth,unmeth) = line
    chrm = chrm.strip()
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    # coord should be in base 0 already, just need it as an int
    coord = int(coord)
    meth = float(meth)
    unmeth = float(unmeth)
    return (chrm,coord,meth,unmeth)


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["aggregateFile=","methBedGraph=","percentageBedGraph=","coverageBedGraph="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    percentageBedGraph = None
    coverageBedGraph = None 
    
    for o, a in opts:
        if o=="--aggregateFile":
            infile = a
            print "Bismark", a
        elif o=="--methBedGraph":
            methBedGraph = a
            print "Meth Bed Graph", a
        elif o=="--percentageBedGraph":
            percentageBedGraph = a
            print "Percentage Bed Graph", a
        elif o=="--coverageBedGraph":
            coverageBedGraph = a
            print "Coverage Bed Graph", a

    assert infile != None
    assert methBedGraph != None
    assert percentageBedGraph != None
    assert coverageBedGraph != None
    
    starttime = datetime.now()
    
    minReadsInAtLeastOneSampleForPercentage = 10
    
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    
    methBed = csv.writer(open(methBedGraph,"w"),delimiter='\t')
    percentBed = csv.writer(open(percentageBedGraph,"w"),delimiter='\t')
    coverageBed = csv.writer(open(coverageBedGraph,"w"),delimiter='\t')

    
    def methScore(methylated,unmethylated):
        methScore = methylated - unmethylated
        # log the methScore (note you can't log a 0 or a negative number so do some fixing for that
        if methScore != 0:
            if methScore > 0:
                methScore = math.log(methScore)
            else:
                methScore = -1*math.log(abs(methScore))
        return methScore
    

    for line in methfile:
        (chrm,coord,meth,unmeth) = parseMethLine(line)
        
        methBed.writerow([chrm,coord,coord+1,methScore(meth,unmeth)])  

        # only write the percentages if at least minReads number of reads
        if meth + unmeth > minReadsInAtLeastOneSampleForPercentage:
            methPercent = (float(meth)/float(meth+unmeth)) * 100.0
            percentBed.writerow([chrm, coord, coord+1, methPercent])
            
        coverageBed.writerow([chrm,coord,coord+1,meth+unmeth])
                  

    print datetime.now()-starttime                    