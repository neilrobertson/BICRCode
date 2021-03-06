'''
Created on 19 Aug 2011

@author: mcbryan
'''

import getopt
import sys
import csv

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["one=", "two=", "combine="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    oneloc = None
    twoloc = None
    combineloc = None
    
    for o, a in opts:
        if o=="--one":
            oneloc = a
            print "One", a
        elif o=="--two":
            twoloc = a
            print "Two", a
        elif o=="--combine":
            combineloc = a
            print "Combine", a
    
    assert oneloc != None
    assert twoloc != None
    assert combineloc != None
    
    left = csv.reader(open(oneloc, "r"), delimiter="\t")
    right = csv.reader(open(twoloc, "r"), delimiter="\t")
    
    combine = csv.writer(open(combineloc,"w"),delimiter='\t')      
    
    def extractLine(line):
        (chr,coord,methylated,unmethylated) = line
        chr = chr.strip()
        # coord should be in base 0 already, just need it as an int
        coord = int(coord)
        methylated=int(methylated)
        unmethylated=int(unmethylated)
        if not chr.startswith("chr"):
            chr = "chr"+chr
        return (chr,coord,methylated,unmethylated)
    
    
    # get the first data points
    
    lchr,lcoord,lmeth,lunmeth = extractLine(left.next())
    rchr,rcoord,rmeth,runmeth = extractLine(right.next())
    
    leftFinished = False
    rightFinished = False
    
    
    def spoolBoth():
        global lchr,lcoord,lmeth,lunmeth,leftFinished
        global rchr,rcoord,rmeth,runmeth,rightFinished
        
        combine.writerow([lchr, lcoord, lmeth, lunmeth, rmeth,runmeth])
        
        try:    
            lchr,lcoord,lmeth,lunmeth = extractLine(left.next())
        except StopIteration:
            leftFinished = True    
        try:            
            rchr,rcoord,rmeth,runmeth = extractLine(right.next())
        except StopIteration:
            rightFinished = True 
        
    def spoolLeft():
        global lchr,lcoord,lmeth,lunmeth,leftFinished
        
        combine.writerow([lchr, lcoord, lmeth,lunmeth,0,0])
            
        try:    
            lchr,lcoord,lmeth,lunmeth = extractLine(left.next())
        except StopIteration:
            leftFinished = True          
    
    def spoolRight():
        global rchr,rcoord,rmeth,runmeth,rightFinished
        
        combine.writerow([rchr, rcoord, 0,0, rmeth,runmeth])
        
        try:            
            rchr,rcoord,rmeth,runmeth = extractLine(right.next())
        except StopIteration:
            rightFinished = True    
    
    while not leftFinished or not rightFinished:
   
        if leftFinished:
            spoolRight()
        elif rightFinished:
            spoolLeft()
        else:
            if lchr == rchr and lcoord == rcoord:
                # we have a match, print it
                spoolBoth()
            elif lchr < rchr:
                spoolLeft()
            elif lchr > rchr:
                spoolRight()   
            # we are on the same chromosome but we have skipped one (or more)
            elif lcoord < rcoord:
                # we have a coord in left that we havnt reached in right yet
                # we therefore want to deal with all the left stuff until we catch up
                spoolLeft()
            else:
                # we have a coord in right that we havnt reached in left yet
                spoolRight()
