'''
Created on 20 Mar 2014

@author: johncole
'''

import scipy.stats
import numpy
import getopt
import sys
import csv
from itertools import izip_longest
from itertools import izip

def unpackLine(chrm,coord,*data):
    return chrm.strip(), int(coord), [int(d) for d in data]

def extractLine(line):
    chrm,coord,data = unpackLine(*line)
    if not chrm.startswith("chr"):
        chrm = "chr"+chrm
    return (chrm,coord,data)

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def pooledMethCounts(data,covs,compare):
    
    methdata = []
    
    for i in compare:
        for cov, (meth,unmeth) in izip(covs,grouper(2,data)): # extract in pairs of meth/unmeth
            if cov == i:
                methdata.append(meth)
                methdata.append(unmeth)
                
    return methdata




if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","cov=","compare=","window=","out=", "minCoverage="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    windowSize = 2000
    infile = None
    covs = None
    compare = None
    outpath = None
    minimumCoverage = 0
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "in file:", a
        elif o=="--cov":
            covs = a.split(",")
            print covs
        elif o=="--compare":
            compare = a.split(",")
            print compare
        elif o=="--window":
            windowSize = int(a)
            print "window size:", a
        elif o=="--out":
            outpath = a
            print "out file:", a
        elif o=="--minCoverage":
            minimumCoverage = int(a)
            print  "Minimum accepted coverage:", a
                
    assert outpath != None
    assert infile != None
    assert covs != None
    assert compare != None
    
    currentWindowPos = None
    currentChr = None
    dataCount = None
    invalidWindows = 0
    cs = 0
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    outfile = open(outpath,"w")

    #We print a header row:
    headerString = "Chr\tStart\tStop\tCpGs_Tested"
    for i in compare:
        headerString += "\t" + i
    headerString += "\tmean\tstdDev\tp-chi.sq"

    outfile.write(headerString + "\n")
    

    #Method to collect and collate the data for a window and print:
    def writeWindow(chrm,pos,epos):
        global dataCount
        global windows
        global outfile
        global invalidWindows
        global cs
        
        #String to be printed
        outputString = str(chrm) +"\t" + str(pos) + "\t" + str(epos+1) + "\t" + str(cs)
        
        #We get the values for the first treatment:
        windowData = pooledMethCounts(dataCount,covs,compare)
        
        windowPercent = []
        methList = []
        unmethList = []
        
        #Test for an unprintable window, if not calculate % meths
        for x in xrange(0,len(windowData),2):
            for y in xrange(1,len(windowData),2):
                if ((windowData[x] + windowData[y]) == 0 or (windowData[x] + windowData[y]) < minimumCoverage)  and y-x ==1:
                    invalidWindows += 1
                    return
                elif  y-x == 1:
                    windowPercent.append((float(windowData[x])/float(windowData[x]+windowData[y])) * 100.0)
                    outputString += "\t" + str(windowPercent[-1])
                    methList.append(windowData[x])
                    unmethList.append(windowData[y])

        #Calculate the mean,Standard deviation and Chi.Sq:
        arr = numpy.array(windowPercent)
        obs = numpy.array([methList, unmethList]).T
        try:
            chi2,p,dof,expected = scipy.stats.chi2_contingency(obs) 
        except ValueError:
            p = 1.0

        outputString += "\t" + str(numpy.mean(arr,axis=None)) + "\t" + str(numpy.std(arr, axis=0)) + "\t" + str(p)
       
        #We print the row
        outfile.write(outputString + "\n")
        
        dataCount = None
        cs = 0
        

    lastCoordCounted = None

    #Read the file and generate windows, calls window to be printed:
    for line in methfile:
        
        chrm,coord,data = extractLine(line)

        if currentChr == None:
            currentChr = chrm
            currentWindowPos = coord
            lastCoordCounted = coord
        else:
            # always write a window if we change chrm
            if chrm != currentChr:
                writeWindow(currentChr,currentWindowPos,lastCoordCounted)
                currentChr = chrm
                currentWindowPos = coord
                print currentChr
                
            # write a window if the next cpg is outside the current window
            elif coord - currentWindowPos > windowSize:
                writeWindow(currentChr,currentWindowPos,lastCoordCounted)
                currentChr = chrm
                currentWindowPos = coord
                
            # otherwise dont write any window (the current cpg is in the window)
        
        ####
        # add data to existing data count
        ####
        if dataCount == None:
            dataCount = data
        else:
            for index,d in enumerate(data):
                dataCount[index] += d
        lastCoordCounted = coord
        
        cs += 1
        
    # after we write a window there is always at least one unmeasured coord after it so we have one more window to write
    writeWindow(currentChr,currentWindowPos,lastCoordCounted)
    
    print invalidWindows
    outfile.flush()
    outfile.close()