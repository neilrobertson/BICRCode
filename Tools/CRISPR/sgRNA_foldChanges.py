'''
Created on 7 Jan 2016

@author: neilrobertson
'''

import getopt,sys
import numpy as np
import scipy.stats as scistats
from collections import defaultdict



def getColumns(inFile, delim="\t", header=True):
        cols = {}
        indexToName = {}
        
        for lineNum, line in enumerate(inFile):
            if lineNum == 0:
                headings = line.split(delim)
                i = 0
                for heading in headings:
                    heading = heading.strip()
                    if header:
                        cols[heading] = []
                        indexToName[i] = heading
                    else:
                        # in this case the heading is actually just a cell
                        cols[i] = [heading]
                        indexToName[i] = i
                    i += 1
            else:
                cells = line.split(delim)
                i = 0
                for cell in cells:
                    cell = cell.strip()
                    cols[indexToName[i]] += [cell]
                    i += 1 
                    
        return cols, indexToName
    
    
def calculatePerSampleReadTotals(columns, controls, treatments):
    print "Calculating sample read totals..."
    totalReadsDict = {}
    
    for rep in controls + treatments: 
        columns[rep] = map(int, columns[rep])
        totalReadsDict[rep] = sum(columns[rep])
        columns[rep] = np.array(columns[rep]).astype(np.float64)
        
    return totalReadsDict


def buildsgRNAMap(columns):
    print "Building sgRNA map..."
    sgRNAs_perGene = defaultdict(list) 
    gene_sgRNA_map = {} 
    sgRNAOrder = []
    geneOrder = []
      
    for i in range(0, len(columns["sgRNA"])):
        sgRNAOrder.append(columns["sgRNA"][i])
        geneOrder.append(columns["Gene"][i])
        gene_sgRNA_map[columns["sgRNA"][i]] = columns["Gene"][i]
        sgRNAs_perGene[columns["Gene"][i]].append(columns["sgRNA"][i])
        
    seen = set()
    seen_add = seen.add
    geneOrder = [x for x in geneOrder if not (x in seen or seen_add(x))]
        
    return sgRNAOrder, geneOrder, gene_sgRNA_map, sgRNAs_perGene


def buildDataMap(columns):
    dataMap = defaultdict(list)
    print "Building data map..."
    
    for i in range(0, len(columns["sgRNA"])):
        current_sgRNA = columns["sgRNA"][i]
        controlRatios = []
        treatmentRatios = []
        for rep in controls: controlRatios.append(np.divide(columns[rep][i], np.float64(totalReadsDict[rep])))
        for rep in treatments: treatmentRatios.append(np.divide(columns[rep][i], np.float64(totalReadsDict[rep])))
        dataMap[current_sgRNA].append(controlRatios)
        dataMap[current_sgRNA].append(treatmentRatios)
        
    return dataMap



if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["sgRNACounts=", "outputPrefix=", "control=", "treatment="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
countsFilename = None
outputFilename = None   
controls = None
treatments = None 
    
DELIMITER = "\t"
    
for o, a in opts:
    if o=="--sgRNACounts":
        countsFilename = a
    elif o=="--outputPrefix":
        outputFilename = a
    elif o=="--control":
        controls = a
    elif o=="--treatment":
        treatments = a
        
assert countsFilename, "No sgRNA counts filename provided. Exiting"
assert outputFilename, "No output filename provided. Exiting."
assert controls, "No control column names provided. Exiting."
assert treatments, "No treatment column names provided. Exiting."

controls = controls.strip().split(",")
treatments = treatments.strip().split(",")

assert len(controls) > 0
assert len(treatments) > 0
        
with open(countsFilename, "r") as inputFile:
    print "Opened file %s for processing..." % (countsFilename.strip().split(r"/")[-1])
    columns, indexToName = getColumns(inputFile)
    print "File read."
    
orderedHeaders = []
for i in range(0, len(indexToName.keys())): orderedHeaders.append(indexToName[i])
for rep in controls + treatments: assert rep in orderedHeaders, "%s is not contained in data file headers. Exiting." % (rep)
controls, treatments = tuple(controls), tuple(treatments) 
print "File contains provided replicate names and contains %s sgRNAs." % (str(len(columns["sgRNA"])))  


totalReadsDict = calculatePerSampleReadTotals(columns, controls, treatments)
    
sgRNAOrder, geneOrder, gene_sgRNA_map, sgRNAs_perGene = buildsgRNAMap(columns)

dataMap = buildDataMap(columns)
    
with open(outputFilename + ".perRNA.csv", "w") as output:
    print "Processing and writing output per sgRNA..."
    headerLine = ["sgRNA", "Gene", "control_geomMean", "treatment_geomMean", "foldChange"]
    output.write(DELIMITER.join(headerLine) + "\n")
    
    for sgRNA in sgRNAOrder:
        gene = gene_sgRNA_map[sgRNA]
        data = dataMap[sgRNA]

        controlMean = scistats.mstats.gmean(data[0])
        treatmentMean = scistats.mstats.gmean(data[1])
        
        foldChange = np.divide(treatmentMean, controlMean)
        
        outputLine = [sgRNA, gene, str(controlMean), str(treatmentMean), str(foldChange)]
        output.write(DELIMITER.join(outputLine) + "\n")



with open(outputFilename + ".perGene.csv", "w") as output:
    print "Processing and writing output per gene..."
    
    headerLine = ["Gene", "Available_sgRNAs", "Represented_sgRNAs", "Mean_foldChange"]
    output.write(DELIMITER.join(headerLine) + "\n")
    
    for gene in geneOrder:
        
        gene_sgRNAs = sgRNAs_perGene[gene]
        number_sgRNAs = len(gene_sgRNAs)
        print gene, str(number_sgRNAs), gene_sgRNAs
        
        controlData = []
        treatmentData = []
        for sgRNA in gene_sgRNAs:
            #print sgRNA
            data = dataMap[sgRNA]
            controlData.append(data[0])
            treatmentData.append(data[1])
        
        controlExperiments = zip(*controlData)
        available_Experiments =  len(controlExperiments)
        available_sgRNAs = len(controlData)
        
        for i in range(0, len(controlData)):
            contData = controlData[i]
            cData = []
            for x in contData:
                if x > 0:
                    cData.append(x)
            if len(cData) > 0:
                controlData[i] = scistats.mstats.gmean(cData)
            else:
                controlData[i] = 0
        for i in range(0, len(treatmentData)):
            treatData = treatmentData[i]
            tData = []
            for x in treatData:
                if x > 0:
                    tData.append(x)
            if len(tData) > 0:
                treatmentData[i] = scistats.mstats.gmean(tData)
            else:
                treatmentData[i] = 0
        
        foldChanges = []
        print controlData
        print treatmentData

        if len(treatmentData) == len(controlData):
            for i in range(0, len(controlData)):

                treatmentMean = treatmentData[i]
                controlMean = controlData[i]
                foldChange = np.divide(treatmentMean, controlMean)
                foldChanges.append(foldChange)
        print foldChanges
        if len(foldChanges) > 0:
            foldChanges = np.array(foldChanges).astype(np.float64)
            
            for x in np.nditer(foldChanges, op_flags=['readwrite']):
                x[...] = np.log2(x)
            print foldChanges
            
            foldChanges = foldChanges[~np.isinf(foldChanges)]
            foldChanges = foldChanges[~np.isnan(foldChanges)]
            
            represented_sgRNAs = len(foldChanges)
            
            print foldChanges
            
            averageLog2Fold = np.average(foldChanges)
            print averageLog2Fold
            averageFoldChange = np.power(10, averageLog2Fold)
            
            print available_sgRNAs, represented_sgRNAs
            print averageFoldChange
            outputLine = [gene, str(available_sgRNAs), str(represented_sgRNAs), str(averageFoldChange)]
            output.write(DELIMITER.join(outputLine) + "\n")
        
        

    
        
        
        
print totalReadsDict.values()

print columns[controls[0]]

print orderedHeaders

print "COMPLETED!"
        
