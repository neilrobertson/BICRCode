'''
Created on 17 Feb 2015

@author: neilrobertson

This script collapses dyad sites into a single CpG locus in the understanding that functionally the methylation at dyads
sites is similar.

Inputs:
--methFile    The input file you wish to collapse dyads
--diadFile    The output file with collapsed dyads
--genome      The reference genome in fasta format
--format      The file format of the input file.  Currently supports coverage files from bismarkToBedGraph and aggregate files of all sizes.
'''
import sys, csv, getopt, gc, time, os.path
import cPickle as pickle
import sitetobed

CPG_DYAD_SEQ = r"CG"
FORMATS = ("COVERAGE", "AGGREGATE")

'''The outputs of various methylation files vary in how their positions are indexed, for instance bismark outputs with an initial
index of 1, whereas the bedgraph, UCSC and aggregate file formats require an initial index of 0. This in turn requires a a means 
to subtract from the dyad list that is dependent on file format'''
FORMAT_INDEXER = {"COVERAGE":0, "AGGREGATE":1}


def getDyadsInfo(genomeFilename,CPG_DYAD_SEQ):
    genomeFasta = genomeFilename.split("/")[-1]+".pickle"
    print genomeFasta
    if os.path.isfile(genomeFasta):
        print "Serialized dictionary of dyad sites found at location: %s  DESERIALIZING..." % (genomeFasta)
        dyadLoci = pickle.load(open(genomeFasta, "rb"))
        return dyadLoci
    else:
        print "No dyads dictionary found. Creating a new list from: %s" % (genomeFilename)
        dyadLoci = sitetobed.getSitesDict(genomeFilename, CPG_DYAD_SEQ)
        pickle.dump(dyadLoci, open(genomeFasta, "wb"))
        print "Deserializing dictionary..."
        return pickle.load(open(genomeFasta, "rb"))
        

def replicate_block(datalist, n=2):
    for i in xrange(0, len(datalist), n):
        yield datalist[i:i+n]
        
def getMethLocation(line):
    (chrm, pos) = line.split("\t")[:2]
    return chrm, int(pos)

def getMethKey(line, FORMAT_SPECIFIC_POSITION_CHANGE):
    (chrm, pos) = line.split("\t")[:2]
    return chrm+"_"+str(int(pos)+FORMAT_SPECIFIC_POSITION_CHANGE)

def collapseDyad(initialLine, dyadLine, fileFormat):
    if fileFormat == "COVERAGE":
        (chrm1, start1, end1, percentMeth1, meth1, unmeth1) = initialLine.rstrip().split("\t") #@UnusedVariable
        (chrm2, start2, end2, percentMeth2, meth2, unmeth2) = dyadLine.rstrip().split("\t") #@UnusedVariable
        dyadMeth = int(meth1) + int(meth2)
        dyadUnmeth = int(unmeth1) + int(unmeth2)
        dyadPercent = (dyadMeth/(dyadMeth + float(dyadUnmeth)))*100
        return [chrm1, start1, end1, ("%.13f" % dyadPercent), str(dyadMeth), str(dyadUnmeth)] 
    elif fileFormat == "AGGREGATE":
        initialReplicates = list(replicate_block(initialLine.rstrip().split("\t")[2:]))
        dyadReplicates = list(replicate_block(dyadLine.rstrip().split("\t")[2:]))
        outputLine = []
        outputLine.append(initialLine.rstrip().split("\t")[0])
        outputLine.append(initialLine.rstrip().split("\t")[1])
        #may fail if there are different counts of replicates on each line
        for i, rep in enumerate(initialReplicates):
            outputLine.append(str((int(rep[0]) + int(dyadReplicates[i][0]))))
            outputLine.append(str((int(rep[1]) + int(dyadReplicates[i][1]))))
        return outputLine
    else:
        print "Format not supported. Check input file and arguments......."
        return -1
    
def rollbackGtoC(lineParts, fileFormat):
    if fileFormat == "AGGREGATE":
        pos = int(lineParts[1]) - 1
        lineParts[1] = str(pos)
    elif fileFormat == "COVERAGE":
        start = int(lineParts[1]) - 1
        end = int(lineParts[2]) - 1
        lineParts[1] = str(start)
        lineParts[2] = str(end)
    return lineParts

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["methFile=","dyadCorrected=","genome=","format="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

methFilename = None
dyadCorrectedFile = None
genome = None
fileFormat = None

for o, a in opts:
    if o=="--methFile":
        methFilename = a
    elif o=="--dyadCorrected":
        dyadCorrectedFile = a
    elif o=="--genome":
        genome = a
    elif o=="--format":
        fileFormat = a.upper()

assert methFilename != None
assert dyadCorrectedFile != None
assert genome != None
assert fileFormat in FORMATS

print "Obtaining list of CpG dyads from genome: %s" % (genome)      
dyadLoci = getDyadsInfo(genome, CPG_DYAD_SEQ)
print "Total dyads found: %s" % (len(dyadLoci.keys()))


print "*"*25
print "Preparing to open massive file: %s" % (methFilename)
#This is rather memory intensive, however it allows us to manipulate the file easier
openedMethFile = open(methFilename, "r")
methFile = openedMethFile.readlines()
FORMAT_SPECIFIC_POSITION_CHANGE = FORMAT_INDEXER[fileFormat]
print "File opened. Format: %s Format specified change from actual dyad position: %s" % (fileFormat, str(FORMAT_SPECIFIC_POSITION_CHANGE))
outputFile = open(dyadCorrectedFile, "w")
output = csv.writer(outputFile, delimiter="\t")

print "Commencing CpG dyad collapse on %s CpG positions..." % (str(len(methFile)))
combinedCpGs = 0
cytoNoGuanCount = 0
guanNoCytoCount = 0
failed = 0
totalDyads = len(dyadLoci)
prevTime = time.time()
for pos in range(0, (len(methFile)-1)):
    if ((pos+combinedCpGs) % 100000) == 0:
        print "Working on row: %s" % (str(pos))
    if (pos + combinedCpGs + 1) <= (len(methFile)-1):
        methKey = getMethKey(methFile[pos + combinedCpGs], FORMAT_SPECIFIC_POSITION_CHANGE)
        try:
            [dyadChrm, dyadStart, dyadEnd] = dyadLoci[methKey]
            
            methChr, methPosition = getMethLocation(methFile[pos + combinedCpGs])
            nextLineChr, nextLinePos = getMethLocation(methFile[pos + combinedCpGs +1])
            if dyadChrm == methChr and (dyadStart == (methPosition + FORMAT_SPECIFIC_POSITION_CHANGE)) and dyadChrm == nextLineChr and (dyadEnd == (nextLinePos + FORMAT_SPECIFIC_POSITION_CHANGE)):
                #This is a bloody diad
                collapsedLine = collapseDyad(methFile[pos + combinedCpGs], methFile[pos + combinedCpGs + 1], fileFormat)
                output.writerow(collapsedLine)
                combinedCpGs += 1
            else:
                #This corresponds to a location that has been found in the dyad dictionary (ie is a cytosine) but no corresponding G was found.
                output.writerow(methFile[pos + combinedCpGs].rstrip().split("\t"))
                cytoNoGuanCount += 1
        except:
            # This location was not found in the dyad dictionary.  It represents a guanine or a cystosine on the reverse strand.
            # We move this location back to the cytosine on the reference genome.
            cytokey = getMethKey(methFile[pos + combinedCpGs], (FORMAT_SPECIFIC_POSITION_CHANGE - 1))# G to C
            try:
                [dyadChrm, dyadStart, dyadEnd] = dyadLoci[cytokey]
                correctedLine = rollbackGtoC(methFile[pos + combinedCpGs].rstrip().split("\t"), fileFormat)
                output.writerow(correctedLine)
                guanNoCytoCount += 0
            except Exception as inst:
                print inst
                failed += 1
                #print "This is not properly converting the reverse strand cytosines (guanines) to forward cytosines."
    else:
        break

openedMethFile.flush()
openedMethFile.close()
outputFile.flush()
outputFile.close()
gc.collect()

print "Files closed..."
print "CpG DYADS COLLAPSED! In total there were %s dyads collapsed. %s skipped." % (str(combinedCpGs), str(totalDyads - combinedCpGs))
print "Forward cytosine without Guanine: %s" % (str(cytoNoGuanCount))
print "Guanines without initial cytosine stepped back: %s" % (str(guanNoCytoCount))
print "Failed counts: %s" % str(failed)
print "*"*25