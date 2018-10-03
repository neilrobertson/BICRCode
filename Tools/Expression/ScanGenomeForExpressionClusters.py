'''
Created on 13 Feb 2014

Tool to scan a list of genes (with expression fc values) for clusters of nearby genes with a similar expression direction

@author: johncole
'''


import sys
import getopt
from bed.treatment import BedFile
import scipy.stats

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:w:s:c:", ["input=","output=","window=","slide", "chrlengths"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    outfile = None
    chrlengthsFile = None
    window = 1000000
    slide = 100000

    for o, a in opts:
        if o=="-i" or o== "--input":
            infile = a
        if o=="-o" or o=="--output":
            outfile = open(a,"w")
        if o=="-w" or o=="--window":
            window = int(a)
        if o=="-s" or o=="--slide":
            slide = int(a)
        if o=="-c" or o=="--chrlengths":
            chrlengthsFile = open(a).readlines()
    
    assert infile !=None
    assert outfile !=None
    assert chrlengthsFile !=None
    assert isinstance(window, int) 
    assert isinstance(slide, int)
    
    #We load the chromosome sizes into memory
    chrlengths = {}
    for line in chrlengthsFile:
        tokens = line.rstrip().split('\t')
        chrlengths[tokens[0]] = int(tokens[1])
    
    #We determine the number of genes, and the number of genes in each direction
    bedArray = open(infile).readlines()
    
    upregulatedGenes = 0
    downregulatedGenes = 0
    totalGenes = 0
    
    for line in bedArray:
        tokens = line.rstrip().split('\t')
        try:
            if float(tokens[3]) > 0:
                upregulatedGenes +=1
            else:
                downregulatedGenes +=1
            totalGenes +=1
        except ValueError:
            pass
            
    probabilityUp = float(upregulatedGenes)/float(totalGenes)
    probabilityDown = float(downregulatedGenes)/float(totalGenes)
    
    #We load the bed file into memory
    bedTree = BedFile(infile).buildIntervalTree()
    
    #We print a header line
    print >> outfile, "chr" + "\t" + "start" + "\t" + "stop" + "\t" + "NumberOfGenes" + "\t" + "UpregulatedGenes" + "\t" + "DownregulatedGenes" + "\t" + "%Upregulated" + "\t" + "P-Binomial" + "\t" + "ExpressionValues" + "\t" + "EnsemblIDs" 
    
    #We scan the genome
    currentStart = None
    currentChr = None
    currentChrLength = None
    
    for chr in chrlengths:
        
        print chr
        
        currentChr = chr
        currentChrLength = int(chrlengths[chr])
        currentStart = 0
        
        while currentStart < currentChrLength:
            
            expressionValues = ""
            ensemblIDs = ""
            
            upregulatedGenes = 0
            downregulatedGenes = 0
            totalGenes = 0
            
            beds = bedTree.getValuesInRange(chr, currentStart, currentStart + window)
            for bedEntry in beds:
                value = bedEntry["value"]
                ensembl = bedEntry["ensembl"]
                
                #We check that the gene expression is not an NA or inf
                try:
                    if float(value) > 0:
                        upregulatedGenes += 1
                    else:
                        downregulatedGenes += 1
                    totalGenes += 1
                    
                    #We update the output strings
                    expressionValues += value + ","
                    ensemblIDs += ensembl + ","
                except ValueError:
                    pass
            
            #We print only windows with genes:
            if totalGenes >0:
                print >> outfile, chr + "\t" + str(currentStart) + "\t" + str(currentStart + window) + "\t" + str(totalGenes) + "\t" + str(upregulatedGenes)  + "\t" + str(downregulatedGenes) + "\t" + str(float(upregulatedGenes)/float(totalGenes)) + "\t" + str(scipy.stats.binom_test(upregulatedGenes,totalGenes,probabilityUp)) + "\t" + expressionValues + "\t" + ensemblIDs           
                
            currentStart += slide
            
    outfile.close()   
        
            
            