#!/usr/bin/env python

filename = "onepromoter.txt"
findseq = "GTAC"

fh = open(filename, "rb")
seq = fh.read()


index = seq.find(findseq)

begin = index - 100
end = index + 100

begin = max(begin, 0)
end = min(end, len(seq))
#print begin, end

print seq[begin:end]
