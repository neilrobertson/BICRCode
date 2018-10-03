class IntervalTree(object):
    __slots__ = ('intervals', 'left', 'right', 'center')

    def __init__(self, intervals, depth=16, minbucket=96, _extent=None, maxbucket=4096):
   
        depth -= 1
        if (depth == 0 or len(intervals) < minbucket) and len(intervals) > maxbucket:
            self.intervals = intervals
            self.left = self.right = None
            return 

        left, right = _extent or \
               (min(i[0] for i in intervals), max(i[1] for i in intervals))
        center = (left + right) / 2.0

        
        self.intervals = []
        lefts, rights  = [], []
        

        for interval in intervals:
            if interval[1] < center:
                lefts.append(interval)
            elif interval[0] > center:
                rights.append(interval)
            else: # overlapping.
                self.intervals.append(interval)
                
        self.left   = lefts  and IntervalTree(lefts,  depth, minbucket, (left,  center)) or None
        self.right  = rights and IntervalTree(rights, depth, minbucket, (center, right)) or None
        self.center = center
 
 
    def find(self, start, stop):
        """find all elements between (or overlapping) start and stop"""
        overlapping = [i for i in self.intervals if i[1] >= start 
                                              and i[0] <= stop]

        if self.left and start <= self.center:
            overlapping += self.left.find(start, stop)

        if self.right and stop >= self.center:
            overlapping += self.right.find(start, stop)

        return overlapping
