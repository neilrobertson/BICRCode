#!/usr/bin/env python

import tempfile

(fh, filename) = tempfile.mkstemp()
fh = open(filename, "wb")
print >> fh, "Hello World!"

print filename
