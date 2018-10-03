#!/usr/bin/env python
# -*- coding: utf-8 -*-



import locale
output_encoding = locale.getpreferredencoding()

print output_encoding

txt = u"●"

with open("tester.txt", "wb") as fh:
	print >> fh, txt.encode("utf-8")

