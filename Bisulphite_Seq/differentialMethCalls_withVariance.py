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

def getVariance(obs):
    p = 1.0
    try:
        chi2,p,dof,expected = scipy.stats.chi2_contingency(obs) #@UnusedVariable
    except ValueError:
        p = 1.0
    return p


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["combined=","cov=","compare1=","compare2=","window=","out=", "minCoverage=", "group=","excludeVariance="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    windowSize = 2000
    infile = None
    covs = None
    compare1 = None
    compare2 = None
    outpath = None
    minimumCoverage = 0
    includeVariance = True
    groups = []
    
    for o, a in opts:
        if o=="--combined":
            infile = a
            print "in file:", a
        elif o=="--cov":
            covs = a.split(",")
            print covs
        elif o=="--compare1":
            compare1 = a.split(",")
            print compare1
        elif o=="--compare2":
            compare2 = a.split(",")
            print compare2
        elif o=="--window":
            windowSize = int(a)
            print "window size:", a
        elif o=="--out":
            outpath = a
            print "out file:", a
        elif o=="--minCoverage":
            minimumCoverage = int(a)
            print  "Minimum accepted coverage:", a
        elif o=="--group":
            groups.append(a.split(","))
            print  "Group added:", a
        elif o=="--excludeVariance":
            includeVariance = False
            print "Include variance: %s" % (includeVariance)
            
                
    assert outpath != None
    assert infile != None
    assert covs != None
    assert compare1 != None
    assert compare2 != None
    
    currentWindowPos = None
    currentChr = None
    dataCount = None
    invalidWindows = 0
    cs = 0
    methfile = csv.reader(open(infile, "r"), delimiter="\t")
    out = open(outpath, "w")
    outfile = csv.writer(out, delimiter="\t")

    #We print a header row:
    headerString = ["chr", "start", "stop", "CpGs_Tested"]
    #headerString = "Chr\tStart\tStop\tCpGs_Tested"
    for i in compare1:
        headerString.append(str(i))
        #headerString += "\t" + i
    headerString.append("mean")
    headerString.append("stdDev")
    if includeVariance:
        headerString.append("chi.sq")
    #headerString += "\tmean\tstdDev\tp-chi.sq"
    for i in compare2:
        headerString.append(str(i))
        #headerString += "\t" + i
    headerString.append("mean")
    headerString.append("stdDev")
    if includeVariance:
        headerString.append("chi.sq")
    headerString.append("Cov1_meth")
    headerString.append("Cov1_unmeth")
    headerString.append("Cov2_meth")
    headerString.append("Cov2_unmeth")
    headerString.append("p-Fishers")
    #headerString += "\tmean\tstdDev\tp-chi.sq\tCov1_Meth\tCov1_Unmeth\tCov2_Meth\tCov2_Unmeth\tp-Fishers"
    counter = 0
    for i in groups:
        counter += 1
        headerString.append("Group"+str(counter)+"_meth")
        headerString.append("Group"+str(counter)+"_unmeth")
        #headerString += "\tGroup" + str(counter) + "_Meth" + "\tGroup" + str(counter) + "_UnMeth"
    print headerString
    outfile.writerow(headerString)
    

    #Method to collect and collate the data for a window and print:
    def writeWindow(chrm,pos,epos):
        global dataCount
        global windows
        global outfile
        global invalidWindows
        global cs
        
        #String to be printed
        outputString = [str(chrm), str(pos), str(epos+1), str(cs)]
        
        #Stores the pooled data for determining DMRs:
        compare1Meth = 0
        compare1Unmeth = 0 
        compare2Meth = 0
        compare2Unmeth = 0
        
        #We get the values for the first treatment:
        windowData = pooledMethCounts(dataCount,covs,compare1)
        
        windowPercent = []
        methList = []
        unmethList = []
        
        #Test for an unprintable window, if not calculate % meths
        for x in xrange(0,len(windowData),2):
            for y in xrange(1,len(windowData),2):
                if ((windowData[x] + windowData[y]) == 0 or (windowData[x] + windowData[y]) <= minimumCoverage)  and y-x ==1:
                    invalidWindows += 1
                    return
                elif  y-x == 1:
                    windowPercent.append((float(windowData[x])/float(windowData[x]+windowData[y])) * 100.0)
                    outputString.append(str(windowPercent[-1]))
                    #outputString += "\t" + str(windowPercent[-1])
                    methList.append(windowData[x])
                    unmethList.append(windowData[y])
                    compare1Meth += windowData[x]
                    compare1Unmeth += windowData[y]

        #Calculate the mean,Standard deviation and Chi.Sq:
        arr = numpy.array(windowPercent)
        obs = numpy.array([methList, unmethList]).T
        
        p = getVariance(obs)
        

        outputString.append(str(numpy.mean(arr,axis=None)))
        outputString.append(str(numpy.std(arr, axis=0)))
        if includeVariance:
            outputString.append(str(p))
        #outputString += "\t" + str(numpy.mean(arr,axis=None)) + "\t" + str(numpy.std(arr, axis=0)) + "\t" + str(p)

        #We get the values for the second treatment:
        windowData = pooledMethCounts(dataCount,covs,compare2)
        
        windowPercent = []
        methList = []
        unmethList = []
        
        #Test for an unprintable window, if not calculate % meths
        for x in xrange(0,len(windowData),2):
            for y in xrange(1,len(windowData),2):
                if ((windowData[x] + windowData[y]) == 0 or (windowData[x] + windowData[y]) <= minimumCoverage)  and y-x ==1:
                    invalidWindows += 1
                    return
                elif  y-x == 1:
                    windowPercent.append((float(windowData[x])/float(windowData[x]+windowData[y])) * 100.0)
                    outputString.append(str(windowPercent[-1]))
                    #outputString += "\t" + str(windowPercent[-1])
                    methList.append(windowData[x])
                    unmethList.append(windowData[y])
                    compare2Meth += windowData[x]
                    compare2Unmeth += windowData[y]

        #Calculate the mean,Standard deviation and Chi.Sq:
        arr = numpy.array(windowPercent)
        obs = numpy.array([methList, unmethList]).T
        
        p = getVariance(obs)

        outputString.append(str(numpy.mean(arr,axis=None)))            
        outputString.append(str(numpy.std(arr, axis=0)))
        if includeVariance:
            outputString.append(str(p))
        #outputString += "\t" + str(numpy.mean(arr,axis=None)) + "\t" + str(numpy.std(arr, axis=0)) + "\t" + str(p)
        
        #We compare the two treatments:
        oddsratio, pvalue = scipy.stats.fisher_exact([[compare1Meth, compare1Unmeth], [compare2Meth, compare2Unmeth]])
        outputString.append(str(compare1Meth))
        outputString.append(str(compare1Unmeth))
        outputString.append(str(compare2Meth))
        outputString.append(str(compare2Unmeth))
        outputString.append(str(pvalue))
        #outputString += "\t" + str(compare1Meth) + "\t" + str(compare1Unmeth) + "\t" + str(compare2Meth) + "\t" + str(compare2Unmeth) + "\t" + str(pvalue)
       
        #We get the total % methylation at groups:
        if len(groups) > 0:
            for group in groups:
                
                groupMeth = 0
                groupUnmeth = 0
                windowData = pooledMethCounts(dataCount,covs,group)
                
                for x in xrange(0,len(windowData),2):
                    for y in xrange(1,len(windowData),2):
                        if y-x ==1:
                            groupMeth += windowData[x]
                            groupUnmeth += windowData[y]
                outputString.append(str(groupMeth))
                outputString.append(str(groupUnmeth))
                #outputString += "\t" + str(groupMeth) + "\t" + str(groupUnmeth)
        
        #We print the row
        #print outputString
        outfile.writerow(outputString)

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
    out.flush()
    out.close()