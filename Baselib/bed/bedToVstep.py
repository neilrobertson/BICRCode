#!/usr/bin/env python

import os
import sys
import getopt
import csv
import subprocess

from bp import breakpoint

sys.path.append(sys.path[0])

def do_chrm(infile, outfile, withZeros, vstepWidth):
    chr = ''
    currentPosition = 1
    currentValue = 0
    
    waitingDecrements = {}
    
    row=[]
    
    # get the first position
    while len(row)==0 or row[0].startswith("#"):
        row = infile.next()
        if len(row)==0 or row[0].startswith("#"):
            continue
        chr = row[0]
        
        assert int(row[1]) <= int(row[2]),  "Check for Start after End"
        
        nextIncrementPosition = int(row[1]) # start
        waitingDecrements[int(row[2])] = 1 # stop
        
    # loop until end of file (which we will work out in the middle of this and use a break to escape when we are done)
    while (True):
        while currentPosition + vstepWidth - 1 < nextIncrementPosition:
            
            if withZeros == True or not currentValue == 0:
                outfile.writerow([chr, currentPosition, currentValue])
            
            currentPosition += vstepWidth
            
            # check if higher than a decrement position and if so remove it and decrease current value
            for decrementPosition in waitingDecrements.keys(): # note iterating over a copy of the list as we need to modify while looping
                if currentPosition > decrementPosition:
                    currentValue -= waitingDecrements[decrementPosition]
                    del waitingDecrements[decrementPosition]
        
        # we are past the increment position so need to increase current value
        currentValue += 1
            
        # get next line from the bed file
        try:
            success = False
            while not success:
                row = infile.next()
                if not (len(row)==0 or row[0].startswith("#")):
                    success = True
            
            chr = row[0]
            assert int(row[1]) <= int(row[2]),  "Check for Start after End: "+row[1]+", "+row[2]
            nextIncrementPosition = int(row[1]) # start
            if int(row[2]) in waitingDecrements:
                waitingDecrements[int(row[2])] = waitingDecrements[int(row[2])] + 1# stop
            else:
                waitingDecrements[int(row[2])] = 1# stop
        except StopIteration:
            # we have reached the end of the file but we almost certainly still have some decrements to take care of
            while len(waitingDecrements)>0 and currentPosition <= max(waitingDecrements.keys()):
                # output value
                if withZeros == True or not currentValue == 0:
                    outfile.writerow([chr, currentPosition, currentValue])
                
                # decrements at the end
                for decrementPosition in waitingDecrements.keys(): # note iterating over a copy of the list as we need to modify while looping
                    if currentPosition > decrementPosition:
                        currentValue -= 1
                        waitingDecrements[decrementPosition] = waitingDecrements[decrementPosition] - 1
                        if waitingDecrements[decrementPosition] == 0:
                            del waitingDecrements[decrementPosition]
                
                currentPosition += vstepWidth
            return

def sortAndUniq(bedfile):
    outfile = bedfile + ".sorted"
    cmd = "sort -k 1,1 -k 2,2n %s | uniq > %s" % (bedfile, outfile)
    os.system(cmd)
    os.remove(bedfile)

def sortFile(bedfile):
    sortchr = subprocess.Popen(args=["sort", "-k", "1,1", "-k", "2,2n"],  stdin=open(bedfile),  stdout=open(bedfile+".sorted", "w"))
    sortchr.wait()
    os.remove(bedfile)

def doFile(bedFileLoc, outfileLoc, vstepWidth, withZeros, withDupes, cores):
    from os.path import basename
    basefile = basename(bedFileLoc)
    # split file by chromosome
    splitfile = subprocess.Popen(args=["awk", "-F", "\t", '{close(f);f=$1}{print > "temp-"f"-'+basefile+'"}', bedFileLoc])
    splitfile.wait()

    chrs = []
    files = os.listdir(os.getcwd())
    for bedfile in files:
        if bedfile.startswith("temp-chr"):
            chrs.append(bedfile)
            
    # sort each chr
    from multiprocessing import Pool
    p = Pool(cores)
    if withDupes:
        p.map(sortFile,chrs)
    else:
        p.map(sortAndUniq,chrs)

    outfile = open(outfileLoc, "w")
    vstepfile = csv.writer(outfile, delimiter="\t")
    for bedfile in chrs:
            # bin each chr
            infile = csv.reader(file(bedfile+".sorted"), delimiter="\t")
            do_chrm(infile, vstepfile, withZeros, vstepWidth)
            os.remove(bedfile+".sorted")
    outfile.close()

if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:i:o:0c:d", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
        
    vstepWidth = 200 #default
    withZeros = False
    withDupes = False
    cores = 4

    infileLoc = ""
    outfileLoc = ""
    
    for o, a in opts:
        # need a vstep width
        if o=="-s":
            vstepWidth = int(a)
        elif o=="-i":
            infileLoc = a
        elif o=="-o":
            outfileLoc = a
        elif o=="-0":
            withZeros = True
        elif o=="-c":
            cores = int(a)
        elif o=="-d":
            withDupes = True

    if infileLoc == "" or outfileLoc == "":
        print "missing file argument!"
        sys.exit(2)

    doFile(infileLoc, outfileLoc, vstepWidth, withZeros, withDupes, cores)
