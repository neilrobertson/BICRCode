print "Load file"
inFile = open("/mnt/50tb/privatedata/Douglas/BRAF_Bisulphite_hg19/data/pooledMethylation/CpG_context_Vector-Braf-Aggrigate.txt").readlines() #Dir

totals = []
firstLine = True

print "Loaded file into memory"
for line in inFile:
    tokens = line.rstrip().split("\t") 
    count = 0 
    for i in range(3, len(tokens) + 1):
        
        if firstLine:
            totals.append(0)
        
        totals[count] += int(tokens[i - 1])
        count += 1
    firstLine = False
    
print "Totals length = %s" % str(len(totals))    
print "Read complete"
count = 0
output = []
for i in totals:
    if count%2 == 0:
        output.append(i)
    else:
        output.append(i)
        print output[0], output[1], output[0]+output[1], float(output[0])/float(output[0]+output[1])
        output = []
    count += 1    

               
