#!/usr/bin/env python

import re

fname = "/home/pzs/serdar/tagcassette/code/trunk/config/ucscbaseconfig.txt"

fh = open(fname)

contents = fh.read()

chrom = "chr9"
start = 110517481
end = 110527483

newpos = "position %s:%d-%d" % (chrom, start, end)
print re.sub("position chr\w\w?:\d+-\d+", newpos, contents)

