import csv

def replicate_block(datalist, n=2):
    for i in xrange(0, len(datalist), n):
        yield datalist[i:i+n]

delimiter = "\t"

outputFilenames = {0:"/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/Meth_CpG_SerNeg1/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected.txt",
                   1:"/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/Meth_CpG_SerNeg2/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected.txt",
                   2:"/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/Meth_CpG_SerNeg3/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected.txt",
                   3:"/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/Meth_CpG_SerPos1/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected.txt",
                   4:"/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/Meth_CpG_SerPos2/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected.txt",
                   5:"/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/Meth_CpG_SerPos3/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected.txt"}

inputFilename = "/mnt/50tb/privatedata/non-Adams/Oliver.Maddocks/hg19_RRBS/aggregate_duplicates/CpG_context_aggregate.txt.binom.0.01FDR_CoverageCorrected.txt"

#create output files
outputFiles = {}
for key in  outputFilenames.keys():
    outputcsv = csv.writer(open(outputFilenames[key],"w"),delimiter=delimiter)
    outputFiles[key] = outputcsv

with open(inputFilename, "r") as inputFile:
    for line in inputFile:
        lineParts = line.rstrip().split(delimiter)
        outputLines = {}
        for i in outputFiles.keys():
            outputLines[i] = [lineParts[0], lineParts[1]]
        counter = 0
        for replicateData in replicate_block(lineParts[2:]):
            meth, unmeth = replicateData[0], replicateData[1]
            outputLines[counter].append(meth)
            outputLines[counter].append(unmeth)
            counter += 1
        for key in outputFiles.keys():
            outputFiles[key].writerow(outputLines[key])