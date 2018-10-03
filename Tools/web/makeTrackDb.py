#!/usr/bin/env python

import sys
import getopt
import os
import uuid

if __name__ == "__main__":
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "", ["root=","url=","name=","title=", "tracks=","subfolders","rmdup","uniq","bams"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # add executing directory as part of path
    sys.path.append(sys.path[0])

    root = None
    url = None
    title = None
    tracks = None
    subfolders = False
    removedups = False
    uniqalign = False
    bams = False

    for opt, arg in opts:
        if opt == "--root":
            root = arg
        elif opt == "--url":
            url = arg
        elif opt == "--title":
            title = arg
        elif opt == "--tracks":
            tracks = arg
        elif opt == "--subfolders":
            subfolders = True
        elif opt == "--rmdup":
            removedups = True
        elif opt == "--uniq":
            uniqalign = True
        elif opt == "--bams":
            bams = True
    
    assert root != None
    assert url != None
    assert title != None
    assert tracks != None
    
    url = url.rstrip("/")
    
    tracks = tracks.rstrip("/")
    try:
        os.symlink(root, tracks+"/"+title)
    except OSError:
        None
    
    print "track " + title
    print "shortLabel "+title
    print "longLabel "+title
    print "superTrack on none"
    print "#priority 1"
    
    print
    
    print "track wigs."+title
    print "container multiWig"
    print "configurable on"
    print "shortLabel "+title+" Signal"
    print "longLabel "+title+" Signal"
    print "visibility full"
    print "priority 1"
    print "type bigWig 0 90000"
    print "autoScale on"
    print "alwaysZero on"
    print "aggregate none"
    print "superTrack "+ title + " full"
    
    print
    
    i = 0
    
    if subfolders:
        for currentDir, dirs, files in os.walk(root.rstrip("/")):
            
            for name in files:
                rmdup = ".rmdup" in name
                uniq = ".uniq" in name
                
                if ".bigWig" in name and (rmdup if removedups else not rmdup) and (uniq if uniqalign else not uniq):
                    i+=1
                              
                    lane = currentDir[currentDir.rindex("/")+1:].replace(".","_")

                    shortLabel = lane
                    longLabel = lane
                    
                    print "track " + lane + "_" + str(uuid.uuid1())
                    print "type bigWig"
                    print "shortLabel %s" % shortLabel
                    print "longLabel %s" % longLabel
                    print "parent wigs."+title
                    print "visibility full"
                    print "bigDataUrl "+url+"/"+title+"/"+currentDir[len(root):]+"/"+name
                    print "color 0,0,0"
                    print "priority "+str(i)
                    print
        
        for currentDir, dirs, files in os.walk(root.rstrip("/")):
            
            for name in files:
                rmdup = ".rmdup" in name
                uniq = ".uniq" in name
        
                if name.endswith(".bam") and (rmdup if removedups else not rmdup) and (uniq if uniqalign else not uniq):
                    i+=1
                              
                    lane = currentDir[currentDir.rindex("/")+1:].replace(".","_")

                    shortLabel = lane
                    longLabel = lane
                    
                    print "track bam." + lane + "_" + str(uuid.uuid1())                    
                    print "type bam"
                    print "shortLabel %s" % shortLabel
                    print "longLabel %s" % longLabel                    
                    print "superTrack "+ title + " pack"
                    print "visibility dense"
                    print "bigDataUrl "+url+"/"+title+"/"+currentDir[len(root):]+"/"+name
                    print "color 0,0,0"
                    print "priority "+str(i)
                    print
    else:
        for currentDir, dirs, files in os.walk(root.rstrip("/")+"/bigWig"):
            for name in files:
                
                if ".bigWig" not in name:
                    continue
                
                rmdup = ".rmdup" in name
                uniq = ".uniq" in name
                
                if ".bigWig" in name and (rmdup if removedups else not rmdup) and (uniq if uniqalign else not uniq):
                    i+=1
                    
                    lane = name[:name.index(".fastq")].replace(".","_")
                    
                    longLabel = lane.replace("_"," ")
                    
                    if ".lane" in name:
                        shortLabel = name[:name.index(".lane")].replace("."," ")
                    else:
                        shortLabel = longLabel
                    
                    print "track " + lane + "_" + str(uuid.uuid1())
                    print "type bigWig"
                    print "shortLabel %s" % shortLabel
                    print "longLabel %s" % longLabel
                    print "parent wigs."+title
                    print "visibility full"
                    print "bigDataUrl "+url+"/"+title+"/bigWig/"+name
                    print "color 0,0,0"
                    print "priority 1."+str(i)
                    print
