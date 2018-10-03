'''
Created on 19 Mar 2014

@author: johncole
'''

import getopt
import sys

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["input=","output=","goTerms=","goAnnotations=","geneNames="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    tokens = []
    goterms = {}
    goannotations = {}
    geneNames = {}
    
    for opt, a in opts:
        if opt=="--input":
            inFile = open(a).readlines()
        elif opt=="--output":
            outFile = open(a,'w')
        elif opt=="--goTerms":
            goFile = open(a).readlines()
        elif opt=="--goAnnotations":
            goAnnotationFile = open(a).readlines()
        elif opt=="--geneNames":
            geneNamesFile = open(a).readlines()

    #We load the go terms into memory
    for line in goFile:
        if (line.startswith('id: ')):
            tokens = line.rstrip().split(' ')
            goID = tokens[1]
        elif (line.startswith('name: ')):
            tokens = line.rstrip().split(' ')
            goCategroy = tokens[1]
            goterms[goID] = goCategroy
    
    #We load the go annotations into memory
    for line in goAnnotationFile:
        tokens = line.rstrip().split('\t')
        goID = tokens[4] 
        geneName = tokens[2]
        #geneFullName = tokens[9]
        
        if goannotations.has_key(geneName):
            goannotations[geneName].append(goID)
        else:
            goannotations[geneName] = [goID]
            #geneNames[geneName] = geneFullName 
   
    #We load the gene names into memory
    for line in geneNamesFile:
        tokens = line.rstrip().split(',')
        geneName = tokens[1]
        geneFullName = tokens[2]
        geneNames[geneName] = geneFullName 
        print geneName + '\t' + geneFullName
  

    #We scan the input file:
    for line in inFile:
        tokens = line.rstrip().split('\t')
        geneName = tokens[0]
        geneOntology = ''
        
        if geneNames.has_key(geneName):
            geneOntology += geneNames[geneName] + '\t'
        else:
            geneOntology += 'Unknown Gene Name' + '\t'
        
        if goannotations.has_key(geneName):
            for goID in goannotations[geneName]:
                if goterms.has_key(goID):
                    geneOntology += goterms[goID] + ','
                else:
                    geneOntology += 'Unknown GO ID' + ','
        else:
            geneOntology += 'Unknown Gene ID'
        
        outFile.write(geneOntology + '\t' + line.rstrip() + '\n')    
         
    outFile.flush()
    outFile.close()

