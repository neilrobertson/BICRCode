'''
Created on 19 May 2014

@author: johncole
'''


import getopt
import sys

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["input=","output=","candidates="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    tokens = []
    candidates = {}
    tfs = {}
    tfNames = []
    inFile = None
    candidateFile = None
    outFile = None
    totalGenes = 0
    totalCandidates = 0
    
    for opt, a in opts:
        if opt=="--input":
            inFile = open(a).readlines()
        elif opt=="--output":
            outFile = open(a,'w')
        elif opt=="--candidates":
            candidateFile = open(a).readlines()
    
    #We load load in the list candidate genes:   
    for line in candidateFile:
        tokens = line.rstrip().split('\t')
        candidates[tokens[0]] = True
        totalCandidates +=1


    #We iterate through the gene/tf file:
    first = True
    for line in inFile:
        
        totalGenes += 1
        if totalGenes%1000 == 0:
            print totalGenes
        
        tokens = line.rstrip().split('\t')
        
        #For tf columns:
        for i in range(15,len(tokens)):
            #If its the first line we 
            if first:
                tfNames.append(tokens[i])
                temp = [0,0,0,0]
                tfs[tokens[i]] = temp
            else:
                temp = tfs[tfNames[i-15]]
                temp[0] += int(tokens[i])
                if int(tokens[i])  > 0:
                    temp[2] += 1
                if tokens[3] in candidates:
                    temp[1] += int(tokens[i])
                    if int(tokens[i])  > 0:
                        temp[3] += 1
                tfs[tfNames[i-15]] = temp
                
        first = False

outFile.write("TF_Name" + "\t" + "Total_Count" + "\t" + "Total_Count_%" + "\t" + "Candidate_Count" + "\t" + "Candidate_Count_%" + "\t" + "Expected_Count" + "\t")
outFile.write("Total_Containing" + "\t" + "Total_Containing_%" + "\t" + "Candidate_Containing" + "\t" + "Candidate_Containing_%" + "\t" + "Expected_Containing" + "\n")
print totalGenes,totalCandidates

for i in tfs:
    temp = tfs[i]
    outFile.write(i + "\t" + str(temp[0]) + "\t" + str((float(temp[0])/float(totalGenes))*100) + "\t" + str(temp[1]) + "\t" + str((float(temp[1])/float(totalCandidates))*100) + "\t" + str(((float(temp[0])/float(totalGenes))*float(totalCandidates))) + "\t")
    outFile.write(str(temp[2]) + "\t" + str((float(temp[2])/float(totalGenes))*100) + "\t" + str(temp[3]) + "\t" + str((float(temp[3])/float(totalCandidates))*100) + "\t" + str(((float(temp[2])/float(totalGenes))*float(totalCandidates))) + "\n")

