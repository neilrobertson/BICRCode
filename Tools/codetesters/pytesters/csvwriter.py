#!/usr/bin/env python

import csv

fh = open("tester.csv", "wb")
print >> fh, "header\r"

writer = csv.writer(fh, delimiter=" ", lineterminator="\r\n")

rows = [[1,2],[3,4]]
for row in rows:
	fh.write("chr")
	writer.writerow(row)

fh.close()


