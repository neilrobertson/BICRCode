import numpy

class GeneSlice():
    def __init__(self, sliceid, start, end):
        self.sliceid = sliceid
        self.start = start
        self.end = end
        
    def __repr__(self):
        return str(self.sliceid) +":"+ str(self.start) +"-"+ str(self.end)

def formatDistanceinKB(id, interval):
    return str(abs(id*interval)/1000.0)

# splits a region up into a predetermined number of chunks
def sliceRegion(strand, start, end, numbChunks):
    assert end>start, "End should be greater than start"

    assert numbChunks <= 100, "This doesn't currently support more than 100 regions (as regions are named by %'s)"

    #split into chunks
    length = abs(end-start)
    chunkLength = -(-length / numbChunks) # ceiling division of length / number of chunks -> i.e. it rounds up instead of down, takes advantage of python rounding behaviour for negative numbers
    
    if strand=="+":
        chunks = [int(i) for i in numpy.linspace(start,end,numbChunks,endpoint=False)]
    elif strand=="-":
        chunks = [int(i-chunkLength) for i in numpy.linspace(end,start,numbChunks,endpoint=False)]
    else:
        assert False, "Unknown strand direction"
    
    assert len(chunks)==numbChunks,  "Chunks count incorrect "+str(len(chunks))

    slices = []
    keys=[]
    centreId = 0
    for c in chunks:
        key = str(centreId)+"-"+str(centreId+100/numbChunks)+"%"
        centreId+=100/numbChunks
        slices.append(GeneSlice(key, c, c+chunkLength))
        keys.append(key)
    
    return keys, slices


# slice with fixed intervals extra
def sliceGene(gene, geneSlices, upstreamDistance, upstreamInterval, downstreamDistance, downstreamInterval):
    (chr, strand, start, end) = gene
        
    assert end>start, "End should be greater than start"
    
    if strand=='+':
        upstreamKeys, upstreamSlices = slicePoint(start, upstreamDistance, upstreamInterval, 0, 0, strand)
        downstreamKeys, downstreamSlices = slicePoint(end, 0, 0, downstreamDistance, downstreamInterval, strand)
    elif strand=='-':
        upstreamKeys, upstreamSlices = slicePoint(end, upstreamDistance, upstreamInterval, 0, 0, strand)
        downstreamKeys, downstreamSlices = slicePoint(start, 0, 0, downstreamDistance, downstreamInterval, strand)
        
    geneKeys, geneSlices = sliceRegion(strand,start,end,geneSlices)
    
    return upstreamKeys + geneKeys + downstreamKeys, upstreamSlices + geneSlices + downstreamSlices

# slice with an extra percent
def sliceGenePercent(gene, geneSlices, upstreamPercent, downstreamPercent):
    (chr, strand, start, end) = gene
    
    #split into chunks
    length = abs(end-start)
    
    chunkLength = -(-length / geneSlices) # ceiling division of length / number of chunks -> i.e. it rounds up instead of down, takes advantage of python rounding behaviour for negative numbers
    
    # make upstreamPercent+downstreamPercent times more slices (since it's the number of slices of the gene that was specified and an extra bit of percentage)
    if strand=="+":
        chunks = [int(i) for i in numpy.linspace(start-upstreamPercent*length,end+downstreamPercent*length,geneSlices*(1+upstreamPercent+downstreamPercent),endpoint=False)]
    elif strand=="-":
        chunks = [int(i-chunkLength) for i in numpy.linspace(end+upstreamPercent*length,start-upstreamPercent*length,geneSlices*(1+upstreamPercent+downstreamPercent),endpoint=False)]
    else:
        print "Unknown strand direction"
        exit()
        
    assert len(chunks)==geneSlices*(1+upstreamPercent+downstreamPercent),  "Chunks count incorrect "+str(len(chunks))
        
    centreId = -upstreamPercent*100
    keys=[]
    slices = []
    
    for c in chunks:
        key = str(centreId)+"->"+str(centreId+100/geneSlices)+"%"
        centreId+=100/geneSlices
        keys.append(key)
        
        slices.append(GeneSlice(key, c, c+chunkLength))

    return keys, slices

# slice around a point with fixed intervals
def slicePoint(point, upstreamDistance, upstreamInterval, downstreamDistance, downstreamInterval, strand = '+'):
    
    if upstreamDistance != 0:
        upId = upstreamDistance/upstreamInterval
        if strand=="+":
            upstreamChunks = range(point-upstreamDistance, point, upstreamInterval)
        elif strand=="-":
            upstreamChunks = range(point+upstreamDistance-upstreamInterval, point-upstreamInterval, -upstreamInterval)
        else:
            assert False, "Unknown strand direction"
        assert len(upstreamChunks)==(upstreamDistance/upstreamInterval), "Upstream chunks count incorrect "+str(len(upstreamChunks))
    else:
        upstreamChunks = []
        
    if downstreamDistance != 0:
        downId = 0
        if strand == '+':
            downstreamChunks = range(point, point+downstreamDistance, downstreamInterval)
        elif strand == '-':
            downstreamChunks = range(point-downstreamInterval, point-downstreamDistance-downstreamInterval, -downstreamInterval)
        else:
            assert False, "Unknown strand direction"
        assert len(downstreamChunks)==(downstreamDistance/downstreamInterval),  "Downstream chunks count incorrect "+str(len(downstreamChunks))
    else:
        downstreamChunks = []

    slices = []
    keys=[]
    
    # -1 from the end of the gene slice is to minimise overlaps
    for c in upstreamChunks:
        key = "up"+formatDistanceinKB(upId, upstreamInterval)+"-"+formatDistanceinKB(upId-1, upstreamInterval)+"kb"
        upId-=1
        slices.append(GeneSlice(key, min(c,c+upstreamInterval), max(c,c+upstreamInterval)))
        keys.append(key)
    for c in downstreamChunks:
        key="down"+formatDistanceinKB(downId, downstreamInterval)+"-"+formatDistanceinKB(downId+1, downstreamInterval)+"kb"
        downId+=1
        slices.append(GeneSlice(key, min(c, c+downstreamInterval), max(c, c+downstreamInterval)))
        keys.append(key)
    
    return keys, slices