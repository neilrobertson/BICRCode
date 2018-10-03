'''
Created on 9 Apr 2014

@author: johncole
'''

import networkx as nx
import matplotlib.pyplot as plt

respondents = {}
interactors = {}
interactionColumns = {}
interactions = []

G=nx.Graph()

inFile = open("/mnt/50tb/privatedata/Peter/GARNet/Survey Results - Raw.csv").readlines()

#We get the respondents:
counter = 0
for line in inFile:
    if counter > 1:
        tokens = line.rstrip().split(',')
        respondents[tokens[1]] = tokens[0]
    counter += 1

#We get the interactors who are not respondents:
tokens = inFile[1].rstrip().split(',')
counter = 0
for i in tokens:
    if counter > 2:
        interactors[i.split('\t')[0]] = "NaN"
        if i.split('\t')[0] not in respondents:
            respondents[i.split('\t')[0]] = "NaN"
    counter += 1

#We add the nodes:
for key in interactors:
    G.add_node(key)
        
#We at each interaction column we get the interactor:
tokens = inFile[1].rstrip().split(',')
counter = 0
for i in tokens:
    for j in respondents:
        if j in i:
            interactionColumns[counter] = j
    counter += 1

#We get the interactions
lineNumber = 0
for line in inFile:
    tokens = line.rstrip().split(',')
    column = 0
    for interaction in tokens:
        if (interaction != "" and column > 2 and lineNumber > 1):
            #print tokens[1], interactionColumns[column],interaction
            G.add_edge(tokens[1],interactionColumns[column]) 
            G[tokens[1]][interactionColumns[column]]['color'] = 'blue'
        column += 1
    lineNumber += 1

nx.draw_spring(G,iterations=5000)
plt.show()
