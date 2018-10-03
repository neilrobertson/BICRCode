#!/usr/bin/env python

import cStringIO

cs = cStringIO.StringIO()
cs2 = cStringIO.StringIO()

for i in range(255):
	cs.write(chr(i))

for i in range(255):
	cs2.write(chr(i))

cs.write(cs2.getvalue())

print cs.getvalue()
