'''
Created on 11 Feb 2011

@author: mcbryan
'''
import sys
import getopt

unique = True

if __name__ == '__main__':


    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:b:", [])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    for o, a in opts:
        if o=="-a":
            aFileLoc = a
        if o=="-b":
            bFileLoc = a
            
    bFile = open(bFileLoc,"r")        
    b = set()
    
    for line in bFile:
        b.add(line.strip())
    
    aFile = open(aFileLoc,"r")
    for line in aFile:
        if line.strip() in b:
            print line.strip()