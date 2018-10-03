'''
Created on 10 Feb 2014

Aligns entries in a FastA file using BLAT, outputs bed, extended bed and PSL format.

@author: cole
'''

import sys
import getopt
import subprocess
import math
from fasta.FastA import FastAFile

###################
# BLATaligner needs:
# /mnt/50tb/repository/3rdparty/UCSC/gfServer -canStop -stepSize=5 start localhost 6000 /mnt/50tb/publicdata/hg19/assembly/twoBit/hg19.2bit
###################

splitString = lambda v, l: [v[i*l:(i+1)*l] for i in range(int(math.ceil(len(v)/float(l))))]

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:o:", ["fastA=","outputFile="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    outputBed = None
    outputExtBed = None
    outputPSL = None

    for o, a in opts:
        if o=="-f" or o== "--fastA":
            infile = a
        if o=="-o" or o=="--outputFile":
            outputBed = open(a + ".bed","w")
            outputExtBed = open(a + "-extended.bed","w")
            outputPSL = open(a + ".psl","w")
        
    assert infile != None
    assert outputBed != None

    def getBlatLocation(line):
        result = line.split("\t")
        return int(result[0]),result[8],result[13],result[15],result[16],result[18],result[20]
    def getBlatTStarts(line):
        result = line.split("\t")
        return [int(y) for y in result[20].split(",")[:-1]]
    def getBlatQStarts(line):
        result = line.split("\t")
        return [int(y) for y in result[19].split(",")[:-1]]         
    def getBlatBlockSizes(line):
        result = line.split("\t")
        return [int(y) for y in result[18].split(",")[:-1]]          
    
    for entry in FastAFile(infile):
             
        sequence_handle = open("sequence", "w")
        sequence_handle.write(str(entry))
        sequence_handle.close()

        blat = subprocess.Popen(args=["/mnt/50tb/repository/3rdparty/UCSC/gfClient", "-out=psl", "-nohead", "localhost", "6000", "/", "sequence" , "stdout"],stdout=subprocess.PIPE)
        output, errors = blat.communicate()
        blatresults = str(output).splitlines()

        if len(blatresults) > 0:
            
            #PSL format:
            for blatresult in blatresults:
                print >> outputPSL, blatresult
            
            index = 1
            
            #Extended bed format:
            for blatresult in  sorted(blatresults, key=lambda blatresult: -1*int(blatresult.split("\t")[0])):
                matches,strand,chrom,start,stop,blockSizes,blockQStarts = getBlatLocation(blatresult)
                
                size = 0
                for i in getBlatBlockSizes(blatresult):
                    size += float(i)
                
                tStarts = ""
                for targetStart in getBlatTStarts(blatresult):
                    tStarts += str(targetStart-int(start)) + ","
                
                queryidentity = float(matches)/float(len(entry.sequence))
                blocksIdentity = float(matches) / size
                 
                print >> outputExtBed, chrom  +'\t' + str(start)  +'\t' + stop  +'\t' + entry.header  +'\t' + "1000"  +'\t' + strand   +'\t' + start  +'\t' + stop  +'\t' + "255,255,255"  +'\t' + str(len(getBlatTStarts(blatresult)))  +'\t' + blockSizes  +'\t' + tStarts  +'\t' +  str(queryidentity) +'\t' + str(blocksIdentity) +'\t' + str(index) +'\t' + chrom  +'-' + str(start)  +'-' + stop +'-' + str(index) + '-' + str(queryidentity) + '-' + str(blocksIdentity) + '----' + entry.header.replace(' ','_')                 
                index +=1
            
            
            #Single Exon Bed format
            for blatresult in  sorted(blatresults, key=lambda blatresult: -1*int(blatresult.split("\t")[0])):
                matches,strand,chrom,start,stop,blockSizes,blockQStarts = getBlatLocation(blatresult)
                identity = float(matches)/float(len(entry.sequence))
                
                size = 0
                for i in getBlatBlockSizes(blatresult):
                    size += float(i)
                
                queryidentity = float(matches)/float(len(entry.sequence))
                blocksIdentity = float(matches) / size
                
                blockStarts = getBlatTStarts(blatresult)
                blockSizes = getBlatBlockSizes(blatresult)
            
                for start, length in zip(blockStarts, blockSizes): 
                    print >> outputBed, chrom +'\t' + str(start) + '\t' + str(start+length) + '\t' +  str(matches) + '\t' + str(queryidentity)  + '\t' + str(blocksIdentity)  + '\t' + entry.header
    
    outputPSL.close()
    outputExtBed.close()
    outputBed.close()
    
    print "Done!"
    