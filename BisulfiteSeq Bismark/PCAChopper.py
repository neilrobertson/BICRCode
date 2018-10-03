import csv

def replicate_block(datalist, n=2):
    for i in xrange(0, len(datalist), n):
        yield datalist[i:i+n]

delimiter = "\t"
threshold = 20

outputFileName = "/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/PCA/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected_PerReplicate-PercentMeth.CoverageThreshold-20.txt"
inputFilename = "/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected.txt"

outputcsv = csv.writer(open(outputFileName,"w"),delimiter=delimiter)
with open(inputFilename, "r") as inputFile:
    counter = 0
    outputcsv.writerow(["","SerNeg1","SerNeg2","SerNeg3","SerPos1","SerPos2","SerPos3"])
    for line in inputFile:
        outputLine = []
        outputLine.append(str(counter))
        counter += 1
        hasData = True # Test if all replicates have data before printing
        for replicateData in replicate_block(line.rstrip().split(delimiter)[2:]):
            meth, unmeth = float(replicateData[0]), float(replicateData[1])
            if (meth + unmeth) > threshold: outputLine.append(str((meth/(meth+unmeth))))
            else: hasData = False; break
            
        if hasData == True: outputcsv.writerow(outputLine)