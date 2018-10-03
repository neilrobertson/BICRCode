#!/usr/bin/env python


import libxml2
import libxslt
import sys

if (len(sys.argv) != 3):
	print "Usage: "+sys.argv[0]+" xsltfile xmlfile"
	sys.exit(1)

styledoc = libxml2.parseFile(sys.argv[1])
style = libxslt.parseStylesheetDoc(styledoc)
doc = libxml2.parseFile(sys.argv[2])
result = style.applyStylesheet(doc, None)
print result
style.freeStylesheet()
doc.freeDoc()
result.freeDoc()
