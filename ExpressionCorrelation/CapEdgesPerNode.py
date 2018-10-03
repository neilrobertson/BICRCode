'''
Created on 29 Jan 2015

@author: johncole
'''

import sys
import getopt

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["infile=","outfile=","maxEdges=","oneWay="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    inFile = None
    outFile = None
    maxCorrelations = 10
    oneWay = False

    for o, a in opts:
        if o== "--infile":
            inFile = open(a).readlines()
            print "File opened: " + a
        if o== "--outfile":
            outFile = open(a,'w')
            print "Outputting to: " + a
        if o== "--maxEdges":
            maxCorrelations = int(a)
            print "Capping the number of edges per node at: " + a
        if o== "--oneWay":
            oneWay = a
            if oneWay == True:
                print "Printing the top n (one way) edges per gene"    
            if oneWay == False:
                print "Printing the top n reciprocal (two way) edges per gene"                 


    print "#### Loading nodes and Edges ####"
    
    #For each gene, we generate parallel lists (stored in gene keyed dictionaries) of r values and the edge:
    genes = {}
    edges = {}
    rejected = {}
    accepted = {}
    
    for line in inFile:
        geneA,geneB,PCC = line.rstrip().split('\t')[:3]
        if geneA in genes:
            temp = genes[geneA]
            temp.append(float(PCC))
            genes[geneA] = temp      
            temp = edges[geneA]
            temp.append(line)
            edges[geneA] = temp
            
        else:
            genes[geneA] = [float(PCC)]
            edges[geneA] = [line]
            
        if geneB in genes:
            temp = genes[geneB]
            temp.append(float(PCC))
            genes[geneB] = temp
            temp = edges[geneB]
            temp.append(line)
            edges[geneB] = temp
        else:
            genes[geneB] = [float(PCC)]
            edges[geneB] = [line]

    print "#### Generating top n edges per node ####"

    #We test for genes with > n edges:
    for gene in genes:
        
        #Gets the top n edges for each gene in the gene list
        a = genes[gene]
        IndexesOfTopEdges = sorted(range(len(a)), key=lambda i: abs(a[i]))[-maxCorrelations:]

        #Prints the top n edges if we only want a one way comparison
        if oneWay:
            for edgeIndex in IndexesOfTopEdges:
                geneA,geneB,PCC =  edges[gene][edgeIndex].rstrip().split('\t')[:3]
                outFile.write(geneA + '\t' + geneB + '\t' + PCC + '\n')
        
        #We determine edges that are in the top n enriched for each gene, and those that arent
        if not oneWay:
            for edgeIndex in range(0,len(a)):
                if edgeIndex not in IndexesOfTopEdges:
                    rejected[edges[gene][edgeIndex]] = 1
                else:
                    accepted[edges[gene][edgeIndex]] = 1
                    
    #Once all genes have been assesed we print only those edges that have been accepted in at least one gene, and not rejected in any gene
    if not oneWay:
        for edge in accepted:
            if edge not in rejected:
                outFile.write(edge)
           
    outFile.flush()
    outFile.close()    
    
    print "#### Finished ####"