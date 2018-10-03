'''
Created on 9 Apr 2014

@author: johncole
'''

import networkx as nx
import matplotlib.pyplot as plt
import sys
import getopt

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["infile=","outfile=","degrees=","candidate=", "maxCorr=","minP=","k=","figSize=","nodeSize=","fontSize=","outDPI=","iterations=","edgeLabels","dataOut"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    inFile = None
    outFile = None
    degrees = 5
    candidate = []
    maxCorrelations = 25
    minPearson = 0.9
    kValue = 0.02
    figSize = 50
    nodeSize = 1000
    fontSize = 11
    outDPI = 1000
    iterations = 1000    
    labels = False
    dataOut = False
    
    
    for o, a in opts:
        if o== "--infile":
            inFile = open(a).readlines()
        if o== "--outfile":
            outFile = a
        if o== "--degrees":
            degrees = int(a)
        if  o== "--candidate":
            candidate.append(a)
        if o== "--maxCorr":
            maxCorrelations = int(a)
        if o== "--minP":
            minPearson = float(a)
        if o== "--k":
            kValue = float(a)
        if o== "--figSize":
            figSize = int(a)
        if o== "--nodeSize":
            nodeSize = int(a)
        if o== "--fontSize":
            fontSize = int(a)
        if o== "--outDPI":
            outDPI = int(a)
        if o== "--iterations":
            iterations = int(a)
        if o== "--edgeLabels":
            labels = True
        if o== "--edgeLabels":
            dataOut = True
    
    
    G=nx.Graph()
    nodes = {}
    edges = {}
    edgeColours = []
    edgeLabels = {}
    edgeList= []
    
    for i in candidate:
        nodes[i] = ""
        print i
    
    #We get stats on the input list:
    genes = {}
    correlations = []
    for line in inFile:
        tokens = line.rstrip().split('\t')
        if tokens[0] in genes and abs(float(tokens[2])) > minPearson:
            genes[tokens[0]] = genes[tokens[0]]+1
        elif abs(float(tokens[2])) > minPearson:
            genes[tokens[0]] = 1
        if tokens[1] in genes and abs(float(tokens[2])) > minPearson:
            genes[tokens[1]] = genes[tokens[1]]+1
        elif abs(float(tokens[2])) > minPearson:
            genes[tokens[1]] = 1
    
    for gene in genes:
        correlations.append(genes[gene])
    
    print len(correlations),len(inFile), len(inFile)/len(correlations)
    
    #We get the edges and nodes (the list grows at each degree):
    for degree in range(0,degrees):
        print degree
        temp = {}
        
        for line in inFile:
            tokens = line.rstrip().split('\t')
            
            if ((tokens[0] in nodes or tokens[1] in nodes) and abs(float(tokens[2])) > minPearson):         
                if (genes[tokens[0]] < maxCorrelations and genes[tokens[1]] < maxCorrelations):  
                    edges[line] = tokens
                    temp[tokens[0]] = ""
                    temp[tokens[1]] = ""
        nodes = temp
    
    #We print some stats on edges and nodes:
    print len(nodes)
    print len(edges)
    
    #We create the network graph:
    for node in nodes:
        G.add_node(node)
    for edge in edges:
        G.add_edge(edges[edge][0],edges[edge][1])
        #We generate lists/dicts for colouring and labelling edges:
        #edgeColours.append(float(edges[edge][2]))
        edgeList.append([edges[edge][0],edges[edge][1]]) 
        edgeLabels[edges[edge][0],edges[edge][1]] = round(float(edges[edge][2]), 3)
        
        if float(edges[edge][2]) > 0:
            edgeColours.append(1.0)
        else:
            edgeColours.append(-1.0)
        
    # position is stored as node attribute data for spring layout (k= distance between nodes)
    pos=nx.spring_layout(G,iterations=iterations,dim=2,k=kValue,scale = 1)
    
    # color by path length from center node
    p=nx.single_source_shortest_path_length(G,candidate[0])
    plt.figure(figsize=(figSize,figSize))
    
    # Draw the figure:
    nx.draw_networkx_edges(G,pos,nodelist=[candidate[0]],alpha=0.4,edge_color=edgeColours,edge_cmap=plt.cm.RdBu_r, edgelist = edgeList, width=1.5) #plt.cm.Reds_r
    nx.draw_networkx_nodes(G,pos,nodelist=p.keys(),
                           node_size=nodeSize,
                           node_color=p.values(),
                           cmap=plt.cm.Reds_r)
    nx.draw_networkx_labels(G,pos, font_size = fontSize)
    if labels:
        nx.draw_networkx_edge_labels(G,pos,edge_labels=edgeLabels)
    
    plt.axis('off')
    plt.savefig(outFile, format='pdf',dpi=outDPI)        
    #plt.show()
    
    if dataOut:
        print "done"
    