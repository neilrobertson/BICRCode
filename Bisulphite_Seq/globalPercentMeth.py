print "Load file"
filename = "/mnt/bam01/data/mouse/Bisulfite-seq/Miller_mm9_WGBS_Rapa-CR_Hepatocytes/data/aggregate/Miller.CB15.Rapa.Rep4/CpG_context_aggregate.DyadCorrected"
inFile = open(filename).readlines() #Dir

totals = []
firstLine = True

totalMeth = 0
totalUnmeth = 0
badLines = 0
print "Loaded file into memory: %s" % (filename)
for line in inFile:
    lineParts = line.rstrip().split("\t") 
    if len(lineParts) == 4:
        totalMeth += int(lineParts[2].strip())
        totalUnmeth += int(lineParts[3].strip())
    else:
        badLines +=1
print "BadLines: " + str(badLines)
print totalMeth, totalUnmeth, totalMeth+totalUnmeth, float(totalMeth)/float(totalMeth+totalUnmeth)