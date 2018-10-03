#!/usr/bin/env python

from subprocess import *
import re

#cmd = "/usr/bin/ipcress /home/pzs/primerdesign/primerdesignpython/trunk/ipcresstmp.txt --sequence /home/pzs/genebuilds/human/*.fasta"
cmd = "/usr/bin/exonerate -q /home/pzs/primerdesign/primerdesignpython/trunk/exoneratetmp.txt -t /home/pzs/genebuilds/human/*.fasta --percent 80"
p = Popen(cmd, shell=True, bufsize=100, stdout=PIPE)
output, dummy = p.communicate()

print "found", len(output.split("\n")), "lines"
