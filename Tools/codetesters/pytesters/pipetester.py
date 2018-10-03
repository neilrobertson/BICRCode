#!/usr/bin/env python

from subprocess import *
import re

chipotle = "/home/pzs/codetesters/perltesters/chipotle_tied_array_version.pl"
infile = "/home/pzs/histone/compositeprofile/trunk/reports/regions/CTCF_reg0.txt"

cmd = "%s --infile %s" % (chipotle, infile)
p = Popen(cmd, shell=True, bufsize=100, stderr=PIPE)
retcode = p.wait()
print retcode
output = p.stderr.read()

reobj = re.search("There were (\d+) peaks", output)
print reobj.group(1)
