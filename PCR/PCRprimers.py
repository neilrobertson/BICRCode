'''
Created on 18 Jan 2011

@author: mcbryan
'''

import shutil
import sys
import getopt
from sequence.genome import Genome
from bed.treatment import BedFile
import subprocess
from Bio.Emboss.Applications import Primer3Commandline
from Bio.Emboss import Primer3
import math
from Bio import SeqIO
from filesystem.mkdir import makeDirectory
from csvfile.indexedcsv import IndexedCSV

build = "hg18"
ucscNames = IndexedCSV("/home/mcbryan/mount/publicdata/hg18/ucsc/ucsc-genes-xref.csv")

###################
# isPcr needs:
# ~/mount/repository/3rdparty/UCSC/gfServer -stepSize=5 start localhost 20000 ~/mount/publicdata/hg18/assembly/twoBit/hg18.2bit
###################

splitString = lambda v, l: [v[i*l:(i+1)*l] for i in range(int(math.ceil(len(v)/float(l))))]

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "r:d:o:", ["regions=","dna=","outputFolder=","productsizes="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    outputfolder = None
    dna = False
    
    # pair of (size, range)
    # size of 0 = no desired product size
    desiredproduct = [
                      (90,"75-130"),
                      (100,"100-200"),
                      (150,"100-200"),
                      (200,"100-200"),
                      ]
    
    for o, a in opts:
        if o=="-r" or o== "--regions":
            infile = a
            regions = BedFile(infile)
        if o=="-d" or o== "--dna":
            infile = a
            regions = SeqIO.parse(open(infile,"r"),"fasta")
            dna = True
        if o=="--productsizes":
            desiredproduct = []
            for size in a.split(","):
                optimal,range = size.split(":")
                desiredproduct.append((int(optimal),range))
            
        if o=="-o" or o=="--outputFolder":
            outputfolder = a
        
    assert infile != None
    assert outputfolder != None
    
    def getBlatLocation(line):
        result = line.split("\t")
        return int(result[0]),result[13],int(result[15]),int(result[16]),result[8]

    def getBlatQStarts(line):
        result = line.split("\t")
        return [int(y) for y in result[19].split(",")[:-1]]
        

    exonboundaries = False
    blockStarts = []
    
    genome = Genome(genomeBuild = build)
    
    makeDirectory(outputfolder)
    
    shutil.copy("arrow-down.gif", outputfolder+"/arrow-down.gif")
    shutil.copy("arrow-none.gif", outputfolder+"/arrow-none.gif")
    shutil.copy("arrow-up.gif", outputfolder+"/arrow-up.gif")
    shutil.copy("sortable.css", outputfolder+"/sortable.css")
    shutil.copy("sortable.js", outputfolder+"/sortable.js")
    
    with open(outputfolder+"/index.html","w") as indexFile:
        
        print >> indexFile, """
            <html><body>
            <head>
            <style type="text/css">
            td{font-size:small;}
            th{font-size:small;}
            </style>
            </head>
        """
        
        regioncount = 0
        
        for region in regions:
            
            regioncount += 1        
            
            regionfilename = str(regioncount)+".html"
            
            with open(outputfolder + "/" + regionfilename,"w") as outputFile:
            
                print >> outputFile, """
                    <html><body>
                    <head>
                    <style type="text/css">
                    td{font-size:small;}
                    th{font-size:small;}
                    </style>
                    <script type="text/javascript" src="sortable.js"></script>
                    </head>
                    <link rel="stylesheet" type="text/css" href="sortable.css" />
                """
                # Sortable JS: http://yoast.com/code/sortable-table/current/sortable.js
                
                print >> outputFile, '<center><table class="index">'
                for productsize,productrange in desiredproduct:
                    print >> outputFile, '<tr><td><a href="#'+str(productsize)+'-'+productrange+'">Jump to product size : '+str(productsize) + ' / '+ productrange + '</a></td></tr>'
                print >> outputFile, '</table></center>'
                
                regioncertain = True
                if dna == False:
                    chr = region.chr
                    start = region.start
                    stop = region.stop
                    
                    if "strand" in regions.keys:
                        strand = region.data["strand"]
                        assert strand in ['+','-']
                    else:
                        strand = '+'
                    
                    sequence = genome.getSequence(chr,(start,stop))
                    
                    if strand == '-':
                        sequence = genome.reverse(genome.complement(sequence))
                    
                    name = chr + ":" + str(start) + "-" + str(stop)
                    print >> outputFile, '<h1><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+build+'&position='+name+'">'+name+"</a></h1>"
                else:
                    name = region.name
                    
                    if name in ucscNames:
                        name = name + " ("+ucscNames[name]["geneSymbol"]+")"
                    
                    sequence = region.seq.tostring()
                    
                    sequence_handle = open("sequence", "w")
                    sequence_handle.write(">"+region.name+"\n"+sequence+"\n")
                    sequence_handle.close()
                    
                    blat = subprocess.Popen(args=["/mnt/50tb/repository/3rdparty/UCSC/gfClient", "-out=psl", "-nohead", "localhost", "20000", "/", "sequence" , "stdout"],stdout=subprocess.PIPE)
                    output, errors = blat.communicate()
                        
                    blatresults = str(output).splitlines()
                    
                    chr = None
                    
                    print >> outputFile, "<h1>"+name+"</h1>"
                    
                    if len(blatresults)<=0:
                        print >> outputFile, '<h2>Unknown genomic location</h2>'
                    else:
                        # just use this for now
                        for blatresult in sorted(blatresults, key=lambda blatresult: -1*int(blatresult.split("\t")[0])):
                            
                            matches,blatchr,blatstart,blatstop,blatstrand = getBlatLocation(blatresult)
                            blatresultname = blatchr + ":" + str(blatstart) + "-" + str(blatstop)                        
    
                            if chr == None:
                                chr = blatchr
                                start = blatstart
                                stop = blatstop
                                strand = blatstrand
                                
                            blockStarts = getBlatQStarts(blatresult)
                            
                            print >> outputFile, '<h2><a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+build+'&position='+blatresultname+'">'+blatresultname+" / " + blatstrand + " strand (Blat score:"+str(matches)+"/"+str(len(sequence))+")</a></h2>"
                    
                    if len(blatresults)>1:
                        regioncertain = False
                        print >> outputFile, '<h3><font color="red">Warning: Blat uncertain about location of dna - using first location for pcr OK checks</font></h3>'
                
                
                
                print "Working at : "+name
                
                print >> indexFile, '<h1><a href="'+regionfilename+'">'+name+'</a></h1>'
                
                print >> outputFile, '<h2>Input sequence / region</h2>'
                print >> outputFile, '<table>'
                
                if dna == True:
                    print >> outputFile, '<tr>'
                    print >> outputFile, "<td>" + sequence + "</td>"
                    print >> outputFile, '</tr>'
                elif len(regions.keys) == 0:
                    print >> outputFile, '<tr>'
                    for column in region.data:
                        print >> outputFile, "<td>" + column + "</td>"
                    print >> outputFile, '</tr>'
                else:
                    print >> outputFile, '<tr>'
                    for key in regions.keys:
                        print >> outputFile, '<th>' + key + '</th>'
                    print >> outputFile, '</tr>'
                    print >> outputFile, '<tr>'
                    for key in regions.keys:
                        print >> outputFile, '<td>' + region.data[key] + '</td>'
                    print >> outputFile, '</tr>'
                print >> outputFile, '</table>'
                
                sequencesections = splitString(sequence,80)
                
                print >> outputFile, '<h2>DNA to feed to primer3</h2><pre>'
                for section in sequencesections:
                    print >> outputFile, section
                print >> outputFile, '</pre>'
                
                sequence_handle = open("sequence", "w")
                sequence_handle.write(">name\n"+sequence+"\n")
                sequence_handle.close()
                
                for productsize,productrange in desiredproduct:
                    
                    id = str(productsize)+'-'+productrange 
                    
                    if productsize == 0:
                        print >> outputFile, '<h2 id="'+id+'">Product range : '+ productrange +'</h2>'
                    else:             
                        print >> outputFile, '<h2 id="'+id+'">Product size : '+ str(productsize)+' (opt) / '+ productrange +'</h2>'
                
                    # eprimer3 --help for details
                    
                    primer3 = Primer3Commandline(cmd="export EMBOSS_PRIMER32_CORE=primer3_core; eprimer32",sequence = "sequence", auto=True, hybridprobe=True)
                    primer3.explainflag = True
                    primer3.mispriminglibraryfile="~/mount/publicdata/"+build+"/primer3/cat_humrep_and_simple.cgi"
                    primer3.numreturn = 50
                    primer3.outfile = "output.pr3"
                    
                    primer3.psizeopt = productsize # optimal primer size, 0 = none
                    primer3.prange = productrange # product size range
                    
                    # primer size
                    primer3.minsize = 18
                    primer3.osize = 20
                    primer3.maxsize = 23
                    
                    # primer temp
                    primer3.mintm = 58
                    primer3.otm = 60
                    primer3.maxtm = 62
                    primer3.maxdifftm = 3 # max temp difference
                    
                    if exonboundaries and len(blockStarts)>1:
                        targets = " ".join([str(i)+","+str(i+1) for i in blockStarts[1:]])
                        primer3.target = '"'+targets+'"'
                    
                    print str(primer3)
                    
                    subprocess.check_call(str(primer3),shell=True)
                    
                    with open("output.pr3","r") as primer3_results:
                    
                        primer3_info = Primer3.read(primer3_results)
                    
                        print >> outputFile, '<table border="1" class="sortable" id="table-'+str(productsize)+'-'+productrange+'"><tr> \
                              <th>isPCR Matches</th> \
                              <th class="unsortable">isPCR Details</th> \
                              <th class="unsortable">UCSC isPCR</th> \
                              <th>Size</th> \
                              <th>Forward_seq</th> \
                              <th>Forward start</th> \
                              <th>Forward length</th> \
                              <th>Forward tm</th> \
                              <th>Forward gc</th> \
                              <th>Reverse seq</th> \
                              <th>Reverse start</th> \
                              <th>Reverse length</th> \
                              <th>Reverse tm</th> \
                              <th>Reverse gc</th> \
                              </tr>'
                        
                        for primer in primer3_info.primers:
                            
                            print >> outputFile, "<tr>"                    
                            
                            # in silico pcr
                            # ~/bin/x86_64/gfPcr localhost 20000 / ATGCCTATGGATGTGTGTGC ACACTGGGAGCAGCAAGATT stdout
                            ispcr = subprocess.Popen(args=["/home/mcbryan/mount/repository/3rdparty/UCSC//gfPcr", "localhost", "20000", "/", primer.forward_seq, primer.reverse_seq, "stdout"],stdout=subprocess.PIPE)
                            output, errors = ispcr.communicate()
                            
                            ispcrresults = str(output).splitlines()
                            matches = []
                            for line in ispcrresults:
                                if line.startswith(">"):
                                    matches.append(line.replace(primer.forward_seq,"").replace(primer.reverse_seq,""))
                            
                            if len(matches)==1:
                                print >> outputFile, '<td><font color="green"><b>1</b></font></td>'
                            else:
                                print >> outputFile, '<td><font color="red"><b>'+str(len(matches))+'</b></font>'
                            
                            
                            print >> outputFile, '<td><table class="pcrdata">'
                            for match in matches:
                                print >> outputFile, '<tr><td>'
                                primerdata, primerlen = match.replace(">","").strip().split(" ")
                                primerchr, primerloc = primerdata.split(":")
                                if "+" in primerloc:
                                    primerstart, primerstop = primerloc.split("+")
                                else:
                                    primerstart, primerstop = primerloc.split("-")
                                primername = primerchr+":"+primerstart+"-"+primerstop
                                print >> outputFile,'<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+build+'&position='+primername+'">'+primerdata+"</a>"
                                if int(primerstart) >= start and int(primerstop) <= stop and int(primerstop)-int(primerstart)+1==primer.size:
                                    print >> outputFile,' <font color="green">OK</font>'
                                else:
                                    print >> outputFile,' <font color="red">X</font>'
                                print >> outputFile, primerlen
                                print >> outputFile,'</td></tr>'                    
                            print >> outputFile, '</table></td>'
                            
                            
                            print >> outputFile, '<td><a href="http://genome.ucsc.edu/cgi-bin/hgPcr?org=Human&db='+build+'&wp_target=genome&wp_f='+primer.forward_seq+'&wp_r='+primer.reverse_seq+'&Submit=submit&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0">PCR</a></td>'
                            
                        
                            # details
                            
                            record = [
                                primer.size,
                                primer.forward_seq,
                                primer.forward_start,
                                primer.forward_length,
                                primer.forward_tm,
                                primer.forward_gc,
                                primer.reverse_seq,
                                primer.reverse_start,
                                primer.reverse_length,
                                primer.reverse_tm,
                                primer.reverse_gc]
                            
                            for entry in record:
                                print >> outputFile, "<td>"+str(entry)+"</td>"
                            
                            print >> outputFile, "</tr>"
                            
                        print >> outputFile, "</table>"
                print >> outputFile, "</body></html>"          
                 
            
        print >> indexFile, "</body></html>"