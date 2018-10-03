'''
Created on 18 Jan 2011

@author: mcbryan
'''

import sys
import getopt
from sequence.genome import Genome
import subprocess
import math
import csv
from filesystem.mkdir import makeDirectory
from csvfile.indexedcsv import IndexedCSV

build = "hg18"
ucscNames = IndexedCSV("/home/mcbryan/mount/publicdata/hg18/ucsc/ucsc-genes-xref.csv")

###################
# isPcr needs:
# ~/bin/x86_64/gfServer -stepSize=5 start localhost 20000 ~/mount/publicdata/hg18/assembly/twoBit/hg18.2bit
###################

splitString = lambda v, l: [v[i*l:(i+1)*l] for i in range(int(math.ceil(len(v)/float(l))))]

if __name__ == '__main__':

    try:
        opts, args = getopt.getopt(sys.argv[1:], "p:o:", ["primers=","outputFolder="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    infile = None
    outputfolder = None
    
    for o, a in opts:
        if o=="-p" or o== "--primers":
            infile = a            
            primerpairs = csv.reader(open(infile, "r"), delimiter='\t')
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
        
    
    genome = Genome(genomeBuild = build)
    
    makeDirectory(outputfolder)
    
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
        
        print >> indexFile, '<table border="1" class="sortable"><tr> \
                              <th>ID</th> \
                              <th>isPCR Matches</th> \
                              <th class="unsortable">isPCR Details</th> \
                              <th class="unsortable">UCSC isPCR</th> \
                              <th>Forward_seq</th> \
                              <th>Forward length</th> \
                              <th>Reverse seq</th> \
                              <th>Reverse length</th> \
                              </tr>'
        
        for id,leftprimer,rightprimer in primerpairs:
            
            id = id.strip()
            leftprimer = leftprimer.strip()
            rightprimer = rightprimer.strip()
            
            #id,leftprimer,rightprimer = primerpairs[primerpair]["id"],primerpairs[primerpair]["leftprimer"].strip(), primerpairs[primerpair]["rightprimer"].strip() 
                         
            print >> indexFile, "<tr><td>"+id+"</td>"
            
            print leftprimer
            print rightprimer
            
            # in silico pcr
            # ~/bin/x86_64/gfPcr localhost 20000 / ATGCCTATGGATGTGTGTGC ACACTGGGAGCAGCAAGATT stdout
            ispcr = subprocess.Popen(args=["/home/mcbryan/bin/x86_64/gfPcr", "localhost", "20000", "/", leftprimer, rightprimer, "stdout"],stdout=subprocess.PIPE)
            output, errors = ispcr.communicate()
            
            ispcrresults = str(output).splitlines()
            
            matches = []
            for line in ispcrresults:
                if line.startswith(">"):
                    matches.append(line.replace(leftprimer.upper(),"").replace(rightprimer.upper(),""))
            
            if len(matches)==1:
                print >> indexFile, '<td><font color="green"><b>1</b></font></td>'
            else:
                print >> indexFile, '<td><font color="red"><b>'+str(len(matches))+'</b></font>'
            
            print >> indexFile, '<td><table class="pcrdata">'
            for match in matches:
                print >> indexFile, '<tr><td>'
                print match
                primerdata, primerlen = match.replace(">","").strip().split(" ")
                primerchr, primerloc = primerdata.split(":")
                if "+" in primerloc:
                    primerstart, primerstop = primerloc.split("+")
                else:
                    primerstart, primerstop = primerloc.split("-")
                primername = primerchr+":"+primerstart+"-"+primerstop
                print >> indexFile,'<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+build+'&position='+primername+'">'+primerdata+"</a>"

                print >> indexFile, primerlen
                print >> indexFile,'</td></tr>'
            print >> indexFile, '</table></td>'
            
            print >> indexFile, '<td><a href="http://genome.ucsc.edu/cgi-bin/hgPcr?org=Human&db='+build+'&wp_target=genome&wp_f='+leftprimer+'&wp_r='+rightprimer+'&Submit=submit&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0">PCR</a></td>'
              
            record = [leftprimer,
                      len(leftprimer)-1,
                      rightprimer,
                      len(rightprimer)-1]
                            
            for entry in record:
                print >> indexFile, "<td>"+str(entry)+"</td>"
              
            print >> indexFile, "</tr>"
            
        print >> indexFile, "</table>"
            
        print >> indexFile, "</body></html>"