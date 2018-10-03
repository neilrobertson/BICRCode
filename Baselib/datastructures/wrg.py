'''
Created on 20 Sep 2010

@author: mcbryan
'''

import bisect
import random
import collections

class WeightedRandomGenerator(object):
    def __init__(self, weights):
        self.keys = []
        self.totals = []
        running_total = 0

        for key in weights:
            running_total += weights[key]
            self.keys.append(key)
            self.totals.append(running_total)

    def next(self):
        rnd = random.random() * self.totals[-1]
        return self.keys[bisect.bisect_right(self.totals, rnd)]

    def __call__(self):
        return self.next()

if __name__ == "__main__":
    chrs = {"chr1":1,"chr5":2,"chr7":3}
    
    wrg = WeightedRandomGenerator(chrs)
    
    count = collections.defaultdict(int)
    for i in range(10000):
        count[wrg()]+=1
    
    print count