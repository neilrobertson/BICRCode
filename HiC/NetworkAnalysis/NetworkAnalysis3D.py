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
from mayavi import mlab
import pymol


kValue = 0.02
weightToken = 6
colourToken = 7
minValue = 0.0
maxValue = 1.0
definedRange = False
delimiter = "\t"

iterations = 500 

def get_HiC_position(line):
    lineParts = line.split(delimiter)
    chr1 = lineParts[0]
    start1 = lineParts[1]
    end1 = lineParts[2]
    return chr1 + ":" + start1 + "-" + end1

def get_HiC_position_contact(line):
    lineParts = line.split(delimiter)
    chr2 = lineParts[3]
    start2 = lineParts[4]
    end2 = lineParts[5]
    return chr2 + ":" + start2 + "-" + end2

def get_HiC_contact_weight(line, weightToken):
    return float(line.split(delimiter)[weightToken])
    
def to_Console(message):
    print "*" * 25
    print message
    print "*" * 25
    
    

inputFilename = "/mnt/50tb/privatedata/Nikolay/HiC_Final/Analysis/Sen_Pro_Hi-C_bedfiles/Networks/enrichment/bedFiles/chrm3/sen-BclI-100kmatTxtlowerTri.bedlowerTri.log.bed.filt.chr3.csv"
imageFilename = ".".join(inputFilename.split(".")[:-1]) + ".Mayavi-Networkx_3D.png"
coOrdsFilename = ".".join(inputFilename.split(".")[:-1]) + ".NetworkX-Co-Ordinates.txt"

to_Console("Reading input file to memory...")
inFile = open(inputFilename).readlines()

H=nx.Graph()

nodeColoursMap = {}
nodePositionMap = {}
positionCounter = 0

to_Console("Adding nodes to visualisation...")
for int, line in enumerate(inFile):
    hicPosition = get_HiC_position(line)
    colourIdentifier = 0
    if len(line.split(delimiter)) > colourToken:
        colourIdentifier = float(line.split(delimiter)[colourToken])
    
    if H.has_node(hicPosition):
        continue
    else:
        H.add_node(hicPosition)
        H.node[hicPosition]["name"] = hicPosition
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


nodeAttributeName = "name"
G=nx.convert_node_labels_to_integers(H,  label_attribute=nodeAttributeName)


to_Console("Generating spring layout...  BOIIING!")
pos=nx.spring_layout(G, iterations=iterations, dim=3, k=kValue, scale = 1, weight="weight")


xyz=np.array([pos[v] for v in sorted(G)])




to_Console("Writing co-ordinates file for node positions...")
with open(coOrdsFilename, "w") as coOrdsFile:
    for coOrd in xyz.tolist(): coOrdsFile.write("\t".join((map(str, coOrd))).strip() + "\n")


exit()
scalars=np.array(G.nodes())+5

mlab.figure(1, bgcolor=(0, 0, 0))
mlab.clf()

to_Console("Drawing plot!")
pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
                    scalars,
                    scale_factor=0.01,
                    scale_mode='none',
                    colormap="Blues", 
                    resolution=20)

pts.mlab_source.dataset.lines = np.array(G.edges())
tube = mlab.pipeline.tube(pts, tube_radius=0.0001)
mlab.pipeline.surface(tube, color=(0.01, 0.01, 0.01))

mlab.savefig(imageFilename)