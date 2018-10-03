import itertools, csv

def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)

with open("/mnt/50tb/privatedata/non-Adams/Brown-Borg.Holly/analysis/RNA-seq/ExpressionProfile/All_DESeq2_Genes_KnownCodingGenes.sig.1.txt", "r") as inputFile:
    
    ensembleDict = {}
    genesList = []
    for line in inputFile:    
        geneID,ensembleID = line.strip().split("\t")
        genesList.append(ensembleID)
        ensembleDict[ensembleID] = geneID
        
    uniqueGenes = list(set(genesList))
    #uniqueGenes = [x for x in genesList if x not in uniqueGenes]
    print len(uniqueGenes)

    with open("/mnt/50tb/privatedata/non-Adams/Brown-Borg.Holly/analysis/RNA-seq/ExpressionProfile/All_Unique_DESeq2_Genes_KnownCodingGenes.sig.1.FPKM.Average.txt", "w") as outputF:
        outputFile = csv.writer(outputF, delimiter="\t")
        header = ["geneName", "DY", "DO", "NY", "NO"]
        outputFile.writerow(header)
        found = []
        with open("/mnt/50tb/privatedata/non-Adams/Brown-Borg.Holly/data/RNA_Seq/FPKM/FPKM_Matrix.csv", "r") as fpkmMatrix:
            for dataLine in fpkmMatrix:
                dataLine = dataLine.strip().split("\t")
                ensembleID = dataLine.pop(0)
                geneName = dataLine.pop(0).upper()
                if ensembleID in uniqueGenes:
                    found.append(ensembleID)
                    outputLine = [ensembleDict[ensembleID]]
                    groupedPoints = grouper(4, dataLine)
                    print dataLine
                    for i in groupedPoints:
                        outputLine.append(str((sum(map(float, i))/4)))
                    print outputLine
                    outputFile.writerow(outputLine)
        
            for i in uniqueGenes:
                if i not in found:
                    print ensembleDict[i]