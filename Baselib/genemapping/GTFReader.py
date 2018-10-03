'''
Created on 25 Jan 2012

@author: mcbryan

https://github.com/chapmanb/bcbb/tree/master/gff

cd gff
sudo python setup.py install


'''

import csv
import os

def parseGTFline(line):
    chrm = line[0]
    source = line[1]
    feature = line[2]
    start = int(line[3]) # note 1 based fully closed
    end = int(line[4]) # note 1 based fully closed
    score = line[5]
    strand = line[6]
    frame = line[7]
    
    attributes = line[8]
    
    return chrm,source,feature,start,end,score,strand,frame,parseAttributes(attributes)

def parseAttributes(unparsedAttributes):
    a = {}
    # split by ; symbol
    # use the csv module so we can take advantage of the quote char rules which makes the code a little simpler
    for attributes in csv.reader([unparsedAttributes],delimiter=";"):
        for attribute in attributes:
            if not attribute.strip()=="":
                for key,value in csv.reader([attribute.strip()],delimiter=" ",quotechar="\""):
                    a[key] = value
    return a


class GTFFile(object):
    def __init__(self,gtffileLoc):
        self.fileLoc = os.path.expanduser(gtffileLoc)
    
    def __iter__(self):
        with open(self.fileLoc, "r") as inputFile:
            rows = csv.reader(inputFile, delimiter="\t")
            
            for row in rows:
                if row[0].startswith("#"):
                    continue # skip header / comments
                yield parseGTFline(row)
                

if __name__ == "__main__":
    
    for line in GTFFile("~/mount/publicdata/hg18/rna-seq/ensembl/Homo_sapiens.NCBI36.54.fixed.gtf"):
        print line
        