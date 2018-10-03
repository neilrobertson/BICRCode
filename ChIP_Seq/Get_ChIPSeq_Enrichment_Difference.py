'''
Created on 1 Dec 2014

@author: johncole
'''

import os
import getopt
import sys
import numpy

#Deals with the run options:
if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["inputList=","summaryOutputFile=","control="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
        
    inFileList = None
    summaryOutputFile = None
    control = False
    
    print "#### Getting Run Options ####"
    
    for opt, o in opts:
        if opt == "--inputList":
            inFileList = open(o).readlines()
        elif opt == "--summaryOutputFile":
            summaryOutputFile = open(o,'w')
        elif opt == "--control":
            control = True
            
    if control:
        summaryOutputFile.write("file\tNo_Control_Diff\tWith_Control_Diff\tNo_Control_Fold\tWith_Control_Fold\tNo_Control_LogFold\tWith_Control_LogFold\n")
        print "Using a control"
    
    print "#### Processing inputs ####"
    
    #Does the work
    for line in inFileList:
        print "Processing inputs: " + line.rstrip()
        inputFile1Path, inputFile2Path, outputFilePath = line.rstrip().split('\t')
        inputFile1 = open(inputFile1Path).readlines()
        inputFile2 = open(inputFile2Path).readlines()
        assert len(inputFile1) == len(inputFile2), "input files of different lengths: " + str(len(inputFile1)) + " vs " + str(len(inputFile2))
        outputFile = open(outputFilePath,'w')
        
        outputFile.write("chr\tstart\tstop\tNo_Control_Diff\tWith_Control_Diff\tNo_Control_Fold\tWith_Control_Fold\tNo_Control_LogFold\tWith_Control_LogFold\n")
        
        totalsList = [0.0,0.0,0.0,0.0,0.0,0.0]
        countsList = [0.0,0.0,0.0,0.0,0.0,0.0]
        
        for i in range(1, len(inputFile1)):
            tokens1 = inputFile1[i].rstrip().split('\t')
            tokens2 = inputFile2[i].rstrip().split('\t')
                        
            #We check that both controls and both ChIPs have data:
            if ( float(tokens1[3]) > 0.0 and float(tokens1[7]) > 0.0 and float(tokens2[3]) > 0.0 and float(tokens2[7]) > 0.0 ):
            
                No_Control_Diff = float(tokens2[6])- float(tokens1[6])
                With_Control_Diff = float(tokens2[14])- float(tokens1[14])
                
                if (float(tokens2[6]) > float(tokens1[6])):
                    No_Control_Fold = (float(tokens2[6]) / float(tokens1[6])) 
                else:
                    No_Control_Fold = -1 / (float(tokens1[6]) / float(tokens2[6]))
                    
                if (float(tokens2[18]) > float(tokens1[18])):
                    With_Control_Fold = (float(tokens2[18]) / float(tokens1[18])) 
                else:
                    With_Control_Fold = -1 / (float(tokens1[18]) / float(tokens2[18]))
                
                No_Control_LogFold = numpy.log2(float(tokens2[6])) - numpy.log2(float(tokens1[6]))
                With_Control_LogFold = numpy.log2(float(tokens2[18])) - numpy.log2(float(tokens1[18]))
                
                outputFile.write(tokens1[0] + "\t" + tokens1[1] + "\t" + tokens1[2] + "\t")
                outputFile.write(str(No_Control_Diff) + "\t")
                outputFile.write(str(With_Control_Diff) + "\t")
                outputFile.write(str(No_Control_Fold) + "\t")
                outputFile.write(str(With_Control_Fold) + "\t")
                outputFile.write(str(No_Control_LogFold) + "\t")
                outputFile.write(str(With_Control_LogFold) + "\t\n")
                
                #Update the totals:
                totalsList[0] += No_Control_Diff
                countsList[0] += 1
                totalsList[1] += With_Control_Diff
                countsList[1] += 1
                totalsList[2] += No_Control_Fold
                countsList[2] += 1
                totalsList[3] += With_Control_Fold
                countsList[3] += 1
                totalsList[4] += No_Control_LogFold
                countsList[4] += 1
                totalsList[5] += With_Control_LogFold
                countsList[5] += 1

            
            else:
                outputFile.write(tokens1[0] + "\t" + tokens1[1] + "\t" + tokens1[2] + "\t" + "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n")
            
        summaryOutputFile.write(str(os.path.basename(outputFilePath)) + "\t" + str(totalsList[0]/countsList[0]) + "\t" + str(totalsList[1]/countsList[1]) + "\t" + str(totalsList[2]/countsList[2]) + "\t" + str(totalsList[3]/countsList[3]) + "\t" + str(totalsList[4]/countsList[4]) + "\t" + str(totalsList[5]/countsList[5]) + "\n")
        
        