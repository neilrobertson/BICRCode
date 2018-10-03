# Author: pzs
# Extended by McBryan to add chrmEnd trimming + individual file handling

import csv
from genemapping.chrmEnds import ChromosomeEnds
import os

def extendOne(fname, extendlen):
    reader = csv.reader(open(fname), delimiter="\t")
    
    writer = csv.writer(open(fname + ".extended", "w"), delimiter="\t")
    for row in reader:
        chr, start, end, strand = row[0], int(row[1]), int(row[2]), row[5]
        
        if strand == "+":
            row[2] = min(ends[chr],start + extendlen)
        elif strand == "-":
            row[1] = max(0,end - extendlen)
        elif strand == ".":
            pass

        writer.writerow(row)

if __name__ == "__main__":
    from glob import glob
    import sys
    if len(sys.argv) != 4:
        print "usage: extendbed.py <beddir> <genomebuild> <extendlength>"
        sys.exit(1)

    path = sys.argv[1]
    genome = sys.argv[2]
    extendlen = int(sys.argv[3])
    assert extendlen >= 0, "Asked to extend by negative amount"
    
    ends = ChromosomeEnds(genome)
    
    if os.path.isdir(path):
        allfiles = glob(path + "/*.bed")
        for fname in allfiles:
            extendOne(fname, extendlen)
    elif os.path.isfile(path):
        extendOne(path, extendlen)