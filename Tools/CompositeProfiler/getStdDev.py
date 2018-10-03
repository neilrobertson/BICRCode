'''
Created on 7 Apr 2014

@author: johncole
'''

import numpy
import sys
import getopt

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err)
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])
    
    infiles = list()
    stdDevTable = []

    
    for opt, a in opts:
        if opt == "-s":
            infiles.append(open(a).readlines())
            print a
      
    for i in range(len(infiles[0])):
        if i%1000 == 0:
            print i
        temp = infiles[0]
        stdDevRow = []
        for j in range(len(temp[0].rstrip().split('\t'))):
            arr = []
            for infile in infiles:
                try:
                    try:
                        arr.append(float(infile[i].rstrip().split('\t')[j]))
                    except IndexError:
                        continue
                except ValueError:
                    continue 
                
            try:
                stdDevRow.append(str(numpy.std(arr, axis=0)))
            except TypeError:
                stdDevRow.append("")
        
        
        stdDevTable.append(stdDevRow)
              
    for i in range(len(temp[0].rstrip().split('\t'))):
        
        arr = []
        total = 0.0
        count = 0
        
        for row in stdDevTable:
            if row[i] != "nan":
                arr.append(row[i])
                total += float(row[i])
                count += 1
        
        if len(arr) > 0:
            print total/count