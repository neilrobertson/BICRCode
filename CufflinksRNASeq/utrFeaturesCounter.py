'''
Created on 21 Mar 2011

@author: mcbryan
'''

import sys
import getopt
from csvfile.genelist import GeneList
from web.RegRNA import UTRelements
import urllib2
from genemapping import Ensembl
from fasta.FastA import FastAFile


if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["ensemblids=", "fasta="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    assembly = "hg19"
    ensemblids = []
    fastaFile = None
    
    for o, a in opts:
        if o=="--ensemblids":
            ensemblids = GeneList(a)
        elif o=="--fasta":
            fastaFile = FastAFile(a)
    
    if len(ensemblids) > 0:
    
        genedata = Ensembl.EnsemblGenes(assembly=assembly)
        
        print "ensemblid", "chr", "start", "stop", "transcripts", "stemloops", "polyas"
        
        for ensemblid in ensemblids:
            
            stemloops = 0
            polyas = 0
            
            for transcriptid in genedata[ensemblid]:
                # http://genome-euro.ucsc.edu/cgi-bin/hgc?db=hg19&g=htcGeneMrna&i=ENST00000314332&o=ensGene&table=ensGene
                
                chrm = genedata[ensemblid][transcriptid].chr
                start = str(genedata[ensemblid][transcriptid].start)
                stop = str(genedata[ensemblid][transcriptid].end)
            
                try:
                    # get predicted mrna from ucsc    
                    url = "http://genome-euro.ucsc.edu/cgi-bin/hgc?db="+assembly+"&g=htcGeneMrna&o=ensGene&table=ensGene&i="+transcriptid+"&c="+chrm+"&l="+start+"&r="+stop
                           
                    req = urllib2.Request(url)
                    response = urllib2.urlopen(req).read()
                    fasta = response[response.index("<TT>")+len("<TT>"):response.index("</TT>")]
    
                    elements = UTRelements(fasta)
                
                    if "Histone 3UTR stem-loop structure (HSL3)" in elements:
                        stemloops+=1
                    if "Polyadenylation Signal (PAS)" in elements:
                        polyas+=1
                except Exception:
                    pass
            
            numbTranscripts = len(genedata[ensemblid])
            
            print ensemblid, chrm, start, stop, numbTranscripts, stemloops, polyas
    

    if fastaFile != None:
        print "header", "stemloops", "polyas"
        
        for fastaRecord in fastaFile.readlines():
            
            stemloops = 0
            polyas = 0
            
            elements = UTRelements(str(fastaRecord))
        
            if "Histone 3UTR stem-loop structure (HSL3)" in elements:
                stemloops=1
            if "Polyadenylation Signal (PAS)" in elements:
                polyas=1

        
            print fastaRecord.header, stemloops, polyas