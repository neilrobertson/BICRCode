import os, csv


##### ENTER BASE ALIGNMENT DIRECTORY #####
baseDir = "/mnt/50tb/privatedata/Annie/Nutlin_CPG_RNA_Seq/Data/aligned"
subDirs = next(os.walk(baseDir))[1]
subDirs = sorted(subDirs)

def getNumbers(string):
    l = []
    for t in string.split():
        try:
            l.append(float(t.strip(')').strip('%').strip('(')))
        except ValueError:
            pass
    return l

with open(baseDir + "/alignmentStats.csv", "w") as outputFile:
    
    outputcsv = csv.writer(outputFile, delimiter="\t")
    
    header = ["", "L-Input", "Mapped", "Mapped%", "Multi-Alignments", "Multi%", "R-Input", "Mapped", "Mapped%", "Multi-Alignments", "Multi%",
              "Overall-Rate", "Pairs-Aligned", "Multi-Alignments", "Multi%", "Discordant-Alignments", "Discordant%", "Concordant-Rate"]

    outputcsv.writerow(header)
    
    for alignDir in subDirs:
        statsFilename = "/".join((baseDir, alignDir, "align_summary.txt"))
        
        statsFileLines = []
        with open(statsFilename, "r") as statsFile:
            line = ""
            for line in statsFile:
                statsFileLines.append(line)
        
            leftInput = [int(s) for s in statsFileLines[1].split() if s.isdigit()][0]  
            
            leftMapped = [int(s) for s in statsFileLines[2].split() if s.isdigit()][0]
            leftMappedPC = getNumbers(statsFileLines[2])[1]
            
            leftMulti = [int(s) for s in statsFileLines[3].split() if s.isdigit()][0]
            leftMultiPC = getNumbers(statsFileLines[3])[1]
            
            rightInput = [int(s) for s in statsFileLines[5].split() if s.isdigit()][0]  
            
            rightMapped = [int(s) for s in statsFileLines[6].split() if s.isdigit()][0]
            rightMappedPC = getNumbers(statsFileLines[6])[1]
            
            rightMulti = [int(s) for s in statsFileLines[7].split() if s.isdigit()][0]
            rightMultiPC = getNumbers(statsFileLines[7])[1]
            
            overallAlignmentPC = getNumbers(statsFileLines[8])[0]
            
            alignedPairs = [int(s) for s in statsFileLines[10].split() if s.isdigit()][0]
            
            multiPair = [int(s) for s in statsFileLines[11].split() if s.isdigit()][0]
            multiPairPC = getNumbers(statsFileLines[11])[1]

            dicordantPair = [int(s) for s in statsFileLines[12].split() if s.isdigit()][0]
            discordantPC = getNumbers(statsFileLines[12])[1]
            
            concordantPC = getNumbers(statsFileLines[13])[0]
            
            output = [alignDir, str(leftInput), str(leftMapped), str(leftMappedPC), str(leftMulti), 
                      str(leftMultiPC), str(rightInput), str(rightMapped), str(rightMappedPC), str(rightMulti), str(rightMultiPC),
                      str(overallAlignmentPC), str(alignedPairs), str(multiPair), str(multiPairPC), str(dicordantPair), str(discordantPC), str(concordantPC)]
    
            outputcsv.writerow(output)
            
