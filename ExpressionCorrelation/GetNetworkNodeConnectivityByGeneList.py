'''
Created on 30 Jan 2015

@author: johncole (Happy birthday John Cole)
'''

'''
### Add method to calculate distance by coordinates...
'''

import networkx as nx
import sys
import getopt
import numpy

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["nodesFile=","edgesFile=","outFile=", "geneLists="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    G=nx.Graph()
    
    nodesFile = None
    edgesFile = None
    outFile = None
    geneListPaths = None
    coordinates = {}
    
    for o, a in opts:
        if o== "--nodesFile":
            nodesFile = open(a).readlines()
            print "Nodes file opened: " + a
        if o== "--edgesFile":
            edgesFile = open(a).readlines()
            print "Edges file opened: " + a
        if o== "--outFile":
            outFile = open(a,'w')
            print "Outputting to: " + a
        if o== "--geneLists":
            geneListPaths = a
            print "Gene list file opened: " + a
    
    
    #We load the nodes and edges
    print "#### Loading nodes and edges ####"
    
    for line in nodesFile:
        node, degree, xCoord, yCoord = line.rstrip().split('\t')
        G.add_node(node)
        
        temp = []
        temp.append(float(xCoord))
        temp.append(float(yCoord))
        coordinates[node] = numpy.array(temp)
        
    print "Nodes loaded"
    
    for line in edgesFile:
        nodeA, nodeB, PCC, pValue = line.rstrip().split('\t')
        G.add_edge(nodeA,nodeB)

    print "Edges Loaded"
        
        
    print "#### Assessing Genes Lists ####"
    
    #We iterate through each gene list:

        
    print "Assesing connectivity of nodes in file:" + str(geneListPaths.rstrip())
    geneList = open(geneListPaths.rstrip()).readlines()
    
    networkGenesInGeneList = []
    distances = []
    geometricDistances = []
    totalDistance = 0
    totalComparisons = 0
    totalGeometricDistance = 0
    
    #We determine which genes in the gene list are found within the network:
    for gene in geneList:
        if G.has_node(gene.rstrip()):
            networkGenesInGeneList.append(gene.rstrip())
    
    if len(networkGenesInGeneList) >= 5 and len(networkGenesInGeneList) <= 500:
    
        #We get the connectivity between each gene pair:
        for i in range(0,len(networkGenesInGeneList)):
            for j in range(i+1,len(networkGenesInGeneList)):
                
                #We get the distance in degress:
                distance = nx.shortest_path_length(G,networkGenesInGeneList[i],networkGenesInGeneList[j])
                totalDistance += distance
                totalComparisons += 1
                distances.append(distance)
                
                #We get the distance in space:
                geneA = coordinates[networkGenesInGeneList[i]]
                geneB = coordinates[networkGenesInGeneList[j]]
                geometricDistance = numpy.linalg.norm(geneA-geneB)
                totalGeometricDistance += geometricDistance
                geometricDistances.append(geometricDistance)
        
        outFile.write(geneListPaths.rstrip() + "\t" + str(len(networkGenesInGeneList)) + "\t" + str(totalComparisons) + "\t" + str(  float(totalDistance) / float(totalComparisons) ) + "\t" + str(numpy.median(numpy.array(distances))) + "\t" + str(  float(totalGeometricDistance) / float(totalComparisons) ) + "\t" + str(numpy.median(numpy.array(geometricDistances))) + "\n")
        
    print "Finished!"
    

