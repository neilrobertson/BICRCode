import os
import collections
from bed.treatment import Bed
from bed.treatment import SimpleBed

class ChromosomeGaps(Bed):
    def __init__(self,  build):
        super(ChromosomeGaps, self).__init__(os.path.expanduser("~/mount/publicdata/"+build+"/gaps."+build))
        
        self.chrmgaps = collections.defaultdict(int)
        self.individualgaps = SimpleBed(os.path.expanduser("~/mount/publicdata/"+build+"/gaps."+build))
        
        for (chr,start,stop) in self.individualgaps:
            self.chrmgaps[chr] += (stop - start)
            
if __name__ == "__main__":
    gaps = ChromosomeGaps("hg18")
    print gaps.chrmgaps["chr1"]
    
    print gaps.chrmgaps["chr2"]