'''
Created on 30 Jan 2015

@author: johncole
'''

import networkx as nx
import sys
import getopt

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["infile=","outPath=","degrees=","seedGene=", "minPCC=","iterations=", "kValue=", "plotAll="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    G=nx.Graph()
    nodes = {}
    edges = {}
    genes = {}
    correlations = []
    noEdges = 0
    maxEdges = 1
    
    inFile = None
    path = None
    outStatsFile = None
    outDegreesFile = None
    outCoordsFile = None
    plotAll = False
    degrees = 5
    minPCC = 0.8
    iterations = 100
    kValue = 0.01
    
    for o, a in opts:
        if o== "--infile":
            inFile = open(a).readlines()
            print "File opened: " + a
        if o== "--outPath":
            path = a
        if o== "--degrees":
            degrees = int(a)
            print "Generating a network plot to " + str(degrees) + " degrees"
        if  o== "--seedGene":
            nodes[a] = 0
            print "Seed gene: " + a + " added"
        if o== "--minPCC":
            minPCC = float(a)
            print "Considering edges with a PCC of at least: " + str(minPCC)
        if o== "--kValue":
            kValue = float(a)
            print "Using a k-vlaue of: " + str(kValue)
        if o== "--iterations":
            iterations = int(a)
            print "Iterating the force directed model " + str(iterations) + " times"
        if  o== "--plotAll":
            plotAll = True
            print "Plotting all nodes"
    
    path = path + "_Deg" + str(degrees) + "_minPCC" + str(minPCC) + "_kValue" + str(kValue) + "_iterations" + str(iterations)
    print "Outputting to: " + path
    outStatsFile = open(path + "_Stats.csv",'w')
    outEdgesFile = open(path + "_Edges.csv",'w')
    outNodesFile = open(path + "_Nodes.csv",'w')
    
    #We load the nodes and edges, and generate some basic stats on the data:
    print "#### Loading nodes and edges ####"
    
    for line in inFile:
        geneA,geneB,PCC = line.rstrip().split('\t')[:3]
        
        #We add/update gene A if it has at least one edge with a PCC within the minPCC
        if geneA in genes and abs(float(PCC)) >= minPCC:
            genes[geneA] = genes[geneA]+1
            noEdges += 1
            
            if genes[geneA] > maxEdges:
                maxEdges = genes[geneA]
            
        elif abs(float(PCC)) >= minPCC:
            genes[geneA] = 1
            noEdges += 1
                        
        #We add/update gene B if it has at least one edge with a PCC within the minPCC
        if geneB in genes and abs(float(PCC)) >= minPCC:
            genes[geneB] = genes[geneB]+1
            
            if genes[geneB] > maxEdges:
                maxEdges = genes[geneB]
                
        elif abs(float(PCC)) >= minPCC:
            genes[geneB] = 1
            
    
    #We print some basic metrics
    print "#### Input Metrics (Filtered) ####"
    print "Number of nodes = " + str(len(genes))
    print "Number of Edges = " + str(noEdges)
    print "Mean edges per node = " + str(float(noEdges) / float(len(genes)))
    print "Maximum edges per node = " + str(maxEdges)
    
    #We output the statistics file:
    print "Outputting statistics to file: " + path + "_Stats.csv"
    for gene in genes:
        outStatsFile.write(gene + "\t" + str(genes[gene]) + "\t\n")   
    outStatsFile.close()   
    
    #We generate the nodes and edges at each degree, and output them into a csv file:
    if not plotAll:
        
        print "#### Obtaining nodes and edges by degree from seed gene(s) ####"
        
        for degree in range(1,degrees +1):
            temp = {}
            
            for line in inFile:
                geneA,geneB,PCC = line.rstrip().split('\t')[:3]
                
                if geneA in nodes and geneB not in nodes and abs(float(PCC)) >= minPCC:      
                    edges[line] = [geneA,geneB,PCC]
                    temp[geneB] = degree
                
                if geneB in nodes and geneA not in nodes and abs(float(PCC)) >= minPCC:      
                    edges[line] = [geneA,geneB,PCC]
                    temp[geneA] = degree
            
            print "Number of nodes in degree " + str(degree) + " = " + str(len(temp))
            nodes = dict(nodes.items() + temp.items())
            
    # We add all the nodes and edges:
    if plotAll:
        for line in inFile:
            geneA,geneB,PCC = line.rstrip().split('\t')[:3]
            
            if abs(float(PCC)) >= minPCC: 
                edges[line] = [geneA,geneB,PCC]
                nodes[geneA] = 0
                nodes[geneB] = 0
            
            
    #We print some stats on edges and nodes:
    print "Total nodes = " + str(len(nodes))
    print "Total edges = " + str(len(edges))
    print "Mean edges per node = " + str( float(len(edges)) / float(len(nodes)))
    
    #We output the edges file:
    print "Outputting edges to file: " + path + "_Edges.csv"
    for edge in edges:
        outEdgesFile.write(edge) 
    outEdgesFile.close()   
    
    #We create the network graph:
    print "#### Generating network graph - boing!!! ####"
    
    for node in nodes:
        G.add_node(node)
    for edge in edges:
        G.add_edge(edges[edge][0],edges[edge][1])
        
    # position is stored as node attribute data for spring layout (k= distance between nodes)
    pos=nx.spring_layout(G,iterations=iterations,dim=2,k=kValue,scale = 1)
    
    print "Outputting nodes and node coordinates to file: " + path + "_Nodes.csv"
    for node in nodes:
        outNodesFile.write(node + "\t" + str(nodes[node]) + "\t" + str(pos[node][0]) + "\t" + str(pos[node][1]) + "\n")
    outNodesFile.close()
    
    print "Finished!"