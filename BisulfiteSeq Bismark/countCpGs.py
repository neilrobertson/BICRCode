'''
Created on 26 Sep 2013

@author: johncole
'''

import re
from sequence.genome import MemoryGenome


#Hardcoded at the moment but we can change this:
motifs = []
regions = []
motifCounts = []

#Opens a list of motifs to scan
motifsFile = open('/mnt/50tb/privatedata/Kirstin/AZA-BS-vs-RNA/TFactors/TF.Methylated.Motifs.csv').readlines()
for line in motifsFile:
    tokens = line.rstrip().split('\t')
    motifs.append(tokens[1])  

#Opens a list of promoter regions
regionsFile = open('/mnt/50tb/privatedata/Kirstin/AZA-BS-vs-RNA/TFactors/Promoter.Regions.bed').readlines()
for line in regionsFile:
    tokens = line.rstrip().split('\t')
    temp = [tokens[0],tokens[1],tokens[2], tokens[3]]
    regions.append(temp) 

#Creates an output
outCounts = open('/mnt/50tb/privatedata/Kirstin/AZA-BS-vs-RNA/TFactors/PromoterMotifs.csv', 'w')
outBed = open('/mnt/50tb/privatedata/Kirstin/AZA-BS-vs-RNA/TFactors/PromoterMotifs.bed', 'w')

#Does the work
mGenome = MemoryGenome("hg18")

for i in regions:
    seq = mGenome.getSequence(i[0],int(i[1]),int(i[2]))
    temp = []
     
    for j in motifs:
        
        #Gets the position of each motif 
        regionMotifs = [(m.start(0), m.end(0)) for m in re.finditer(j, seq.upper())]
        for k in regionMotifs:
            outBed.write(str(i[0]) + "\t" + str(int(i[1])+int(k[0])) + "\t" + str(int(i[1])+int(k[0])+len(j)) + "\t" + j + "\t" + i[3] + "\n")
        
        #Gets the occurences of each motif
        temp.append(len(re.findall(j,seq.upper())))
     
    motifCounts.append(temp)


#Prints the report
outCounts.write('chr\tstart\tstop\tENSEMBL')
for i in range(len(motifs)):
    outCounts.write('\t' + str(motifs[i]))
outCounts.write('\n')
for i in range(len(regions)):
    for j in  range(len(regions[i])):
        outCounts.write(str(regions[i][j]) + "\t")
    for j in  range(len(motifCounts[i])):
        outCounts.write(str(motifCounts[i][j]) + "\t")
    outCounts.write('\n')
    
    
    