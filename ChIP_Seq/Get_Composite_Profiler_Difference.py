'''
Created on 3 Dec 2014

@author: johncole
'''

import getopt
import sys

#Deals with the run options:
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["inputFileList=","output="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
        
        inputFileList = None
    
    print "#### Loading input files #####"
    for opt, o in opts:
        if opt == "--inputFileList":
            inputFileList = open(o).readlines()
        if opt == "--output":
            output = open(o,'w')
            
        
    print "#### Working ####"
    for files in inputFileList:
        
        print "Processing files: " + files.rstrip()
        input1Path,input2Path,name = files.rstrip().split('\t')
        
        input1 = open(input1Path).readlines()
        input2 = open(input2Path).readlines()
        
        #Generates blank lists to store the totals and counts:
        totalsList = [0] * len(input1[0].rstrip().split('\t')[1:])
        countsList = [0] * len(input1[0].rstrip().split('\t')[1:])
        
        for i in range(1,len(input1)):
        
            tokens1 = input1[i].rstrip().split('\t')[1:]
            tokens2 = input2[i].rstrip().split('\t')[1:]
            assert len(tokens1) == len(tokens2)
            
            for j in range(0,len(tokens1)):
                if tokens2[j] != "" and tokens1[j] != "":
                    diff = float(tokens2[j]) - float(tokens1[j])
                    totalsList[j] += diff
                    countsList[j] += 1
                    
        output.write(name)
        for i in range(0,len(totalsList)):
            output.write("\t" + str((totalsList[i]/countsList[i])))
        output.write("\n")      
    
    print "#### Done! ####" 
                
                
            
        
    
        
            
    