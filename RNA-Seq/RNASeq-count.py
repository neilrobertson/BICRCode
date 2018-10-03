import sys
import getopt
import collections
from csvfile.indexedcsv import IndexedCSV
from rpy import *
from genemapping import UCSC
import csv

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["covdesc="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    covdesc = None
    suppliedCovdesc = False
    for o,a in opts:
        if o=="--covdesc":
            # reflects order of arguments,
            # i.e. ./script --covdesc CTRL,CTRL,CTRL,HIRA,HIRA,HIRA s1 s2 s3 s4 s5 s6
            covdesc = a.split(",") 
            suppliedCovdesc = True
    
    if len(args)<2:
        print "Insufficient number of arguments for comparison"
    elif len(args) == 2:
        # need to pool across all replicates and just compare the two we have against each other
        if not suppliedCovdesc:
            covdesc = [args[0],args[1]] # by definition args is two long at this time
    elif len(args) > 2:
        # can pool across types
        assert suppliedCovdesc
    
    
    
    class Count():
        def __init__(self,assigned,min,max):
            self.assigned = assigned
            self.min = min
            self.max = max
        
    genedata = UCSC.UCSCTranscripts(assembly="hg18")
    
    transcript_counts = collections.defaultdict(dict)
    cluster_counts = collections.defaultdict(dict)

    def addCounts(counts,arg,suffix):
        file = IndexedCSV(arg+suffix)
        for row in file:
            counts[row][arg] = Count(file[row]["AssignedCount"],
                                     file[row]["MinCount"],
                                     file[row]["MaxCount"])
    for arg in args:    
        addCounts(transcript_counts,arg,"-transcripts.csv")        
        addCounts(cluster_counts,arg,"-clusters.csv")
    
    def writeSummary(outfilename,counts,args):
        with open(outfilename,'w') as outfile:
            outcsv = csv.writer(outfile,delimiter="\t")
            header = ["id"]
            for arg in args:
                header.append(arg)
            outcsv.writerow(header)
            for key in counts:
                row = [key]
                for arg in args:
                    if arg in counts[key]:
                        row.append(counts[key][arg].assigned)
                    else:
                        row.append(0)
                outcsv.writerow(row)
    writeSummary("transcripts.csv",transcript_counts,args)
    writeSummary("clusters.csv",cluster_counts,args)
    
    # R commands
    
    r.library("DESeq")
    set_default_mode(NO_CONVERSION)
    
    def deseqCounts(inputFile, outputFilePrefix, covdesc,symbolFunction):
        inputData = r.read_table(inputFile,header=True, row_names=1, sep="\t")
        
        cds = r.newCountDataSet(r.round(inputData),covdesc)
        cds = r.estimateSizeFactors(cds)    
        if len(args)==2: 
            cds = r.estimateVarianceFunctions(cds,pool=True)
        else:
            cds = r.estimateVarianceFunctions(cds)
        
        factors = []
        for factor in covdesc:
            if factor not in factors:
                factors.append(factor)
        
        for a in range(len(factors)):
            for b in range(a+1,len(factors)):
                assert a!=b
                res = r.nbinomTest(cds,factors[a],factors[b])
                
                outputFile = outputFilePrefix+factors[a]+"-vs-"+factors[b]+".csv"
                r.write_table(res, file=outputFile,sep="\t",quote=False)
                annotateTable(outputFile,symbolFunction)
    
    def annotateTable(filename,symbolFunction):        
        with open(filename) as inputFile:
            headerSeen = False
            csvFile = csv.reader(inputFile,delimiter="\t")
            with open(filename+"-annotated.csv","w") as outputFile:
                outcsvFile = csv.writer(outputFile,delimiter="\t")
                for row in csvFile:            
                    if not headerSeen:
                        headerSeen = True
                        row.insert(0,"")
                        row.append("Symbol")
                    else:
                        # id in second column
                        row.append(symbolFunction(row[1]))
                        # fix "down" fold changes to be -1/fc
                        if row[5]=="Inf" or row[5]=="-Inf":
                            pass
                        elif row[5]=="0":
                            row[5]="-Inf"
                        else:
                            if row[5]!="NA" and float(row[5])<1:
                                row[5]=-1.0/float(row[5])

                    outcsvFile.writerow(row)

    def transcriptSymbol(transcriptid):
        return str(genedata[transcriptid].geneSymbol)

    def clusterSymbol(symbolid):
        return str(genedata.clusters[symbolid])
                
    
    deseqCounts("transcripts.csv","transcripts.results.", covdesc, transcriptSymbol)
    deseqCounts("clusters.csv","clusters.results.",covdesc, clusterSymbol)
    
    
    # now annotate
    
    
    print "Done"