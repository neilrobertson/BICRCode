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
        opts, args = getopt.getopt(sys.argv[1:],"", ["nodesFile=","edgesFile=","outNodes=","outEdges=", "geneList="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    G=nx.Graph()
    
    nodesFile = None
    edgesFile = None
    outNodesFile = None
    outEgesFile = None
    geneList = None
    outGenes = {}
    outNodes = {}
    coordinates = {}
    
    for o, a in opts:
        if o== "--nodesFile":
            nodesFile = open(a).readlines()
            print "Nodes file opened: " + a
        if o== "--edgesFile":
            edgesFile = open(a).readlines()
            print "Edges file opened: " + a
        if o== "--outNodes":
            outNodesFile = open(a,'w')
            print "Outputting nodes to: " + a
        if o== "--outEdges":
            outEdgesFile = open(a,'w')            
            print "Outputting edges to: " + a
        if o== "--geneList":
            geneList = open(a).readlines()
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
        
        
    print "#### Assessing Connectivity of Genes in  List ####"
    
    #We iterate through each gene list:
    for geneA in geneList:
        
        geneA = geneA.rstrip()
        print "Assesing gene: " + geneA
        for geneB in geneList:
            
            geneB = geneB.rstrip()
            
            try:
                if geneA in G and geneB in G:
                    path = nx.shortest_path(G,geneA,geneB)
                    
                    for node in path:
                        outNodes[node] = coordinates[node]
            except nx.exception.NetworkXNoPath:
                continue
    
    nodes = outNodes.keys()
    for line in edgesFile:
        nodeA, nodeB, PCC, pValue = line.rstrip().split('\t')
        if nodeA in nodes and nodeB in nodes:
            outEdgesFile.write(line)
    
    for line in nodesFile:
        node, degree, xCoord, yCoord = line.rstrip().split('\t')
        if node in nodes:
            outNodesFile.write(line)
               
    print "Finished!"
    

