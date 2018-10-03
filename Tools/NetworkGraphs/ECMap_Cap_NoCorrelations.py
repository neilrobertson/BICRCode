'''
Created on 10 Apr 2014

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
        if o== "--outfile":
            outFile = open(a,'w')
        if o== "--maxEdges":
            maxCorrelations = int(a) 
        if o== "--oneWay":
            oneWay = a         

    #We get parallel lists (stored in dictionaries) of r values and edges:
    genes = {}
    edges = {}
    rejected = {}
    accepted = {}
    
    for line in inFile:
        tokens = line.rstrip().split('\t')
        if tokens[0] in genes:
            temp = genes[tokens[0]]
            temp.append(float(tokens[2]))
            genes[tokens[0]] = temp      
            temp = edges[tokens[0]]
            temp.append(line)
            edges[tokens[0]] = temp
            
        else:
            genes[tokens[0]] = [float(tokens[2])]
            edges[tokens[0]] = [line]
            
        if tokens[1] in genes:
            temp = genes[tokens[1]]
            temp.append(float(tokens[2]))
            genes[tokens[1]] = temp
            temp = edges[tokens[1]]
            temp.append(line)
            edges[tokens[1]] = temp
        else:
            genes[tokens[1]] = [float(tokens[2])]
            edges[tokens[1]] = [line]
    
    #We test for genes with > n edges:
    for gene in genes:
        
        a = genes[gene]
        sortedList = sorted(range(len(a)), key=lambda i: abs(a[i]))[-maxCorrelations:]
        
        if oneWay:
            for i in sortedList:
                tokens =  edges[gene][i].rstrip().split('\t')
                outFile.write(tokens[0] + '\t' + tokens[1] + '\t' + tokens[2] + '\n')
                
        if not oneWay:
            for i in range(0,len(a)):
                if i not in sortedList:
                    rejected[edges[gene][i]] = 1
                else:
                    accepted[edges[gene][i]] = 1
                    
    if not oneWay:
        for edge in accepted:
            if edge not in rejected:
                outFile.write(edge)
           
    outFile.flush()
    outFile.close()    
    
    