'''
Created on 5 Jan 2015

@author: neilrobertson
'''


import networkx as nx
import matplotlib.pyplot as plt
import sys
import getopt


if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["infile=","outfile=","weights=","colours=","min=","max=","dimensions=","k"])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    #Some hard coded graphical parameters:
    figSize = 50
    nodeSize = 50
    fontSize = 11
    outDPI = 1000
    iterations = 500    

    #Parameters
    inFile = None
    outFile = None
    kValue = 0.02
    weightToken = 6
    colourToken = 7
    minValue = 0.0
    maxValue = 1.0
    definedRange = False
    dimensions = "2d"
    
    for o, a in opts:
        if o== "--infile":
            inFile = open(a).readlines()
        if o== "--outfile":
            outFile = a
        if o== "--weights":
            weightToken = int(a)
        if o== "--colours":
            colourToken = int(a)  
        if o== "--min":
            minValue = float(a)   
            definedRange = True
        if o== "--max":
            maxValue = float(a)       
            definedRange = True
        if o== "--dimensions":
            dimensions = a          
        if o== "--k":
            kValue = float(a)                  
           
    assert dimensions in ("2d", "3d")
         
    G=nx.Graph()
    nodeColoursMap = {}
    nodePositionMap = {}
    positionCounter = 0

    print "#### Updating nodes and edges ####"
    
    for line in inFile:
        tokens = line.rstrip().split('\t')
       
        #We append the node:
        if G.has_node(tokens[0] + ":" + tokens[1] + "-" + tokens[2]):
            continue
        else:
            G.add_node(tokens[0] + ":" + tokens[1] + "-" + tokens[2])
            nodeColoursMap[tokens[0] + ":" + tokens[1] + "-" + tokens[2]] = float(tokens[colourToken])
            nodePositionMap[tokens[0] + ":" + tokens[1] + "-" + tokens[2]] = [positionCounter,positionCounter]
            positionCounter +=1
    
    for line in inFile:
        tokens = line.rstrip().split('\t')

        #We append the edge:
        G.add_edge(tokens[0] + ":" + tokens[1] + "-" + tokens[2],tokens[3] + ":" + tokens[4] + "-" + tokens[5], weight = float(tokens[weightToken]))
    
    print "Number of nodes: " + str(G.number_of_nodes())
    print "Number of edges: " + str(G.number_of_edges())  
    
    print "#### Converting colour map to colours ####"
    values = [nodeColoursMap.get(node) for node in G.nodes()]
    
    # position is stored as node attribute data for spring layout (k= distance between nodes)
    print "#### Generating spring layout  BOING! ####"
    pos=nx.spring_layout(G,iterations=iterations,dim=2,k=kValue,scale = 1, weight="weight") #, pos = nodePositionMap
    
    print "#### Drawing nodes ####"
    if definedRange:
        nx.draw_networkx_nodes(G,pos,node_size=nodeSize,node_color=values, cmap="bwr", vmax = maxValue, vmin = minValue)
    else:
        nx.draw_networkx_nodes(G,pos,node_size=nodeSize,node_color=values, cmap="Reds", vmax = max(values), vmin = min(values))
    
    print "#### Drawing edges ####"
    nx.draw_networkx_edges(G,pos,alpha=0.4, width=1.5)


    print "#### Showing figure  Ooh-la-la! ####"
    plt.show()
    
    
    print "#### Printing figure ####"
    #plt.figure(figsize=(figSize,figSize))
    #plt.axis('off')
    #plt.savefig(outFile, format='pdf',dpi=outDPI)        
