'''
Created on 30 Jan 2015

@author: johncole (Happy birthday John Cole)
'''

import networkx as nx
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["nodesFile=","edgesFile=","outPath=","colourByPath=","colourByValue=", "colourEdges=","drawLabels="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    nodeSize = 100
    nodeFontSize = 12
    nodeFont_weight = "normal"
    edgeAlpha = 1
    edgeWidth = 0.5
    
    G=nx.Graph()
    nodes = {}
    edges = {}
    coordinates = {}
    colours = {}
    seedGene = None
    colourEdges = False
    drawLabels = False
    
    nodesFile = None
    edgesFile = None
    colourFile = None
    path = None
    
    for o, a in opts:
        if o== "--nodesFile":
            nodesFile = open(a).readlines()
            print "Nodes file opened: " + a
        if o== "--edgesFile":
            edgesFile = open(a).readlines()
            print "Edges file opened: " + a
        if o== "--outPath":
            path = a
        if  o== "--colourByPath":
            seedGene = a
            print "Seed gene: " + a + " added"
        if  o== "--colourByValue":
            colourFile = open(a).readlines()
            print "Colours file opened: " + a
        if  o== "--colourEdges":
            colourEdges = True
            print "Colouring edges by PCC"
        if  o== "--drawLabels":
            drawLabels = True
            print "Drawing node labels"
                
    path = path
    print "Outputting to: " + path
    
    #We load the nodes and edges
    print "#### Loading nodes and edges ####"
    
    for line in nodesFile:
        node, degree, xCoord, yCoord = line.rstrip().split('\t')
        nodes[node] = int(degree)
        temp = []
        temp.append(float(xCoord))
        temp.append(float(yCoord))
        coordinates[node] = np.array(temp)
        
    print "Nodes loaded"
    
    for line in edgesFile:
        nodeA, nodeB, PCC, pValue = line.rstrip().split('\t')
        edges[line] = [nodeA,nodeB,PCC]

    print "Edges Loaded"
    
    print "#### Drawing Network ####"
    
    #We generate the image:
    for node in nodes:
        G.add_node(node)
    for edge in edges:
        G.add_edge(edges[edge][0],edges[edge][1])
    
    #We plot the nodes:
    if not seedGene == None:
        p=nx.single_source_shortest_path_length(G,seedGene)
        nx.draw_networkx_nodes(G,coordinates,node_size=nodeSize,nodelist=p.keys(),node_color=p.values(),cmap=plt.cm.RdBu_r)
    
    elif not colourFile == None:
        for line in colourFile:
            node, value = line.rstrip().split('\t')
            colours[node] = float(value)
        nx.draw_networkx_nodes(G,coordinates,node_size=nodeSize,nodelist=colours.keys(),node_color=colours.values(),cmap="bwr")  #"bwr","Paired"
         
    else:
        nx.draw_networkx_nodes(G,coordinates,node_size=nodeSize)
    
    #We plot the edges:
    if colourEdges:
        edgeColours = []
        edgeList = []
        for edge in edges:
            edgeColours.append(float(edges[edge][2]))
            edgeList.append([edges[edge][0],edges[edge][1]])
        nx.draw_networkx_edges(G,coordinates,alpha=edgeAlpha, width=edgeWidth, edge_color=edgeColours,edgelist = edgeList,edge_cmap=plt.cm.RdBu_r, edge_vmin = -1.0,edge_vmax = 1.0)
    
    else:
        nx.draw_networkx_edges(G,coordinates,alpha=edgeAlpha, width=edgeWidth)
    
    #We plot the node labels:
    if drawLabels:
        nx.draw_networkx_labels(G,coordinates, font_size = nodeFontSize,font_weight = nodeFont_weight)
    
    plt.show()
    
    print "Finished!"

