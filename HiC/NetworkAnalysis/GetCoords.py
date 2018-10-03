#############################################################################
#
# Hi-C Network Analysis plotting script
#
# Dependencies: numpy, scipy, matplotlib, networkx and most importantly mayavi2 and it's dependencies
#
# Note: I haven't actually managed to get this working in my IDE.  Have been running this python script in the
# shell using the format with input files hard coded:
# $mayavi2 /home/neilrobertson/workspace/HiC/NetworkAnalysis/NetworkAnalysis3D.py
#
# Author: Neil Robertson
#############################################################################

import networkx as nx
import numpy as np
import sys,getopt


kValue = 0.02
weightToken = 6
colourToken = 7
minValue = 0.0
maxValue = 1.0
definedRange = False

iterations = 500 

DELIMITER = "\t"
NODE_ATTRIBUTE_NAME = "name"

def get_HiC_position(line):
    lineParts = line.split(DELIMITER)
    chr1 = lineParts[0]
    start1 = lineParts[1]
    end1 = lineParts[2]
    return chr1 + ":" + start1 + "-" + end1

def get_HiC_position_contact(line):
    lineParts = line.split(DELIMITER)
    chr2 = lineParts[3]
    start2 = lineParts[4]
    end2 = lineParts[5]
    return chr2 + ":" + start2 + "-" + end2

def get_HiC_contact_weight(line, weightToken):
    return float(line.split(DELIMITER)[weightToken])
    
def to_Console(message):
    print "*" * 25
    print message
    print "*" * 25
    
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["input=","output="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
   
inFilename = None
outFilename = None

for opt, a in opts:
    if opt=="--input":
        inFile = str(a)
    elif opt=="--output":
        outFile = str(a)
    
assert inFilename
assert outFilename
    

to_Console("Reading input file to memory...")
inFile = open(inFilename).readlines()

H=nx.Graph()

nodeColoursMap = {}
nodePositionMap = {}
positionCounter = 0

to_Console("Adding nodes to visualisation...")
for line in inFile:
    hicPosition = get_HiC_position(line)
    colourIdentifier = 0
    if len(line.split(DELIMITER)) > colourToken:
        colourIdentifier = float(line.split(DELIMITER)[colourToken])
    
    if H.has_node(hicPosition):
        continue
    else:
        H.add_node(hicPosition)
        H.node[hicPosition][NODE_ATTRIBUTE_NAME] = hicPosition
        nodeColoursMap[hicPosition] = colourIdentifier
        nodePositionMap[hicPosition] = [positionCounter,positionCounter]
        positionCounter +=1
        

to_Console("Adding edges and claves to nodes...")
for line in inFile:
    hicPosition = get_HiC_position(line)
    positionContact  = get_HiC_position_contact(line)
    hicContactWeight = get_HiC_contact_weight(line, weightToken)
    #We append the edge:
    H.add_edge(hicPosition,positionContact, weight = hicContactWeight)
    
    
print "Number of nodes: " + str(H.number_of_nodes())
print "Number of edges: " + str(H.number_of_edges())

colourValues = [nodeColoursMap.get(node) for node in H.nodes()]

G=nx.convert_node_labels_to_integers(H,  label_attribute=NODE_ATTRIBUTE_NAME)

to_Console("Generating spring layout...  BOIIING!")
pos=nx.spring_layout(G, iterations=iterations, dim=3, k=kValue, scale = 1, weight="weight")
xyz=np.array([pos[v] for v in sorted(G)])


to_Console("Writing co-ordinates file for node positions...")
with open(outFilename, "w") as coOrdsFile:
    for x in sorted(G):
        position = pos[x]
        name = G.node[x]["name"] 
        coOrdsFile.write(name + DELIMITER + DELIMITER.join(map(str, position)).strip() + "\n")
