
from bed.treatment import ExtendedBed
from sam.SamFormat import SAMFile
import os
import sys


lads = ExtendedBed(os.path.expanduser("~/mount/privatedata/non-Adams/donahue.greg/LADs.bed"))

alignments = SAMFile(os.path.expanduser(sys.argv[1]))

numbInLAD = 0
numbNotInLAD = 0

currentSeq = None

def previousKey(key,isin,isnotin):
    print key,isin,isnotin,isin/float(isin+isnotin)

for samEntry in alignments:
    
    if samEntry.chrm == "*":
        continue
    
    if currentSeq != None and samEntry.key != currentSeq:
        # print the stats on the previous key
        previousKey(currentSeq,numbInLAD,numbNotInLAD)
        numbInLAD = 0
        numbNotInLAD = 0
    
    currentSeq = samEntry.key
    
    inLAD = len(lads.getValuesInRange(samEntry.chrm, samEntry.start, samEntry.start+len(samEntry.sequence))) > 0
    
    if inLAD:
        numbInLAD += 1
    else:
        numbNotInLAD +=1

# last one
previousKey(currentSeq,numbInLAD,numbNotInLAD)

