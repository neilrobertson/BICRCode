'''
Created on 23 Aug 2010

@author: mcbryan
'''
import sys
import getopt
from sequence.genome import Genome, UnknownChromosomeException
from genemapping.Ensembl import EnsemblGenes
from csvfile.indexedcsv import IndexedCSV
import collections
import csv

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["matrixdir=","jaspar=","tfd=","gene-expression-file="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    
    # we need the matrix dir first
    matrixdir = None
    for o, a in opts:
        if o=="--matrixdir":
            matrixdir = a
            print "Matrix Dir:",matrixdir            
    assert matrixdir!=None
    
    # now load
    matrices = None
    for o,a in opts:
        if o=="--jaspar":
            print "JASPAR:", a
            
            import JASPAR
            j = JASPAR.Jaspar()
            matrices = j.load(a,matrixdir)
            
        elif o=="--tfd":
            print "TFD:",a
            
            import TFD
            t = TFD.TFD()
            matrices = t.load(a,matrixdir)
    assert matrices != None
    
    exprfile = None
    ensemblidcol = "ensemblid"
    upstreamPromotor = 5000
    downstreamPromotor = 1000
    for o,a in opts:
        if o=="--gene-expression-file":
            exprfile = a
    assert exprfile != None
    
    print len(matrices)
    
    # WARNING: everything we just read in are pfm (position frequency), we might need position weight matrices
    # We have choice of JASPER (downloaded), JASPER (with motility) or TFD or combining them
    
    genedata = EnsemblGenes(assembly="hg18")
    genome = Genome(genomeBuild = "hg18")
    
    exprCSV = IndexedCSV(exprfile,keyPos=1)
    
    motifnames = []
    for matrix in matrices:
        motifnames.append(matrix)

    output = csv.writer(open(exprfile+".transcriptionfactors","w"),delimiter='\t')
        
    # header
    header = exprCSV.keys[:] # copy keys
    header.extend(motifnames)
    output.writerow(header)
    
    for testid in exprCSV:
        
        ensembl = exprCSV[testid][ensemblidcol]
            
        if ensembl not in genedata:
            # RNA-Seq is not a single unique gene, skip
            continue
        
        try:
            promotorstart, promotorstop = genedata[ensembl].getPromotorRegion(upstreamPadding = upstreamPromotor, downstreamPadding = downstreamPromotor)
            
            sequence = genome.getSequence(genedata[ensembl].chr,promotorstart,promotorstop).upper()
            
            motifcount = collections.defaultdict(int)
            
            # we loop over motifnames to make sure it always comes out in the same order as the header
            for motifname in motifnames:
                foundmotifs = matrices[motifname].find(sequence,matrices[motifname].max_score()*0.9)
                for start,end,strand,seq in foundmotifs:
                    if "CG" in seq:
                        motifcount[motifname] += 1
            
            datarow = exprCSV.rows[testid][:] # copy row from expression file
            # we loop over motifnames to make sure it always comes out in the same order as the header
            for motifname in motifnames:
                datarow.append(str(motifcount[motifname]))
            output.writerow(datarow)
            
            print testid
            
        except UnknownChromosomeException:
            continue
