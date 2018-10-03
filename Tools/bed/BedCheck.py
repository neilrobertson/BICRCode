'''
Tool to check that the positions in a given bed file are within the confines of a given genome.
Accepts a bed file, a list of chromosome lengths, and an output path
'''

if __name__ == '__main__':
    
    import getopt
    import sys
    import csv
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["out=","in=","chrms=", "filter="])
    except getopt.GetoptError, err:
        print str(err)
        sys.exit(2)
    
    bedFile = None
    outFile = None
    chrmsFile = None
    filterFlag = False
    
    for o, a in opts:
        if o=="--out":
            outFile = open(a,'w') 
        elif o=="--in":
            inFile = open(a).readlines()
        elif o=="--chrms":
            chrmsFile = open(a).readlines()
        elif o=="--filter":
            filterFlag = True
       
    chrms = {}
        
    for line in chrmsFile:
        chr,length = line.rstrip().split('\t')
        chrms[chr] = length
    
    for line in inFile:
        tokens = line.rstrip().split('\t')
        
        if not tokens[0].startswith("chr"):
                tokens[0] = "chr" + tokens[0]
                        
        if(float(tokens[1]) < 0):
            tokens[1] = str(0)
            print tokens[1]
            
        if(chrms.has_key(tokens[0])):
            if(float(tokens[2]) > float(chrms[tokens[0]])):
                tokens[2] = chrms[tokens[0]]
            
            outFile.write(tokens[0])
            for i in range(1,len(tokens)):
                outFile.write('\t'+ tokens[i])
            outFile.write('\n')   
        
        elif (filterFlag == False):
            outFile.write(tokens[0])
            for i in range(1,len(tokens)):
                outFile.write('\t'+ tokens[i])  
            outFile.write('\n')    
        
    