#!/usr/bin/env python

import urllib
import re
import os
import sys

# gets a screenshot from UCSC of a specifid bed file
# a lot of the work is done in the bed file itself
# and a session file, which tells UCSC what tracks should display and how.
# Upgrades of this code might include exposing a few of those session parameters
# through the object

baseurl = "http://genome.ucsc.edu/cgi-bin/"
pagename = "hgTracks"
convertcmd = "/usr/bin/convert"


class UCSCTools:
	def __init__(self, params=None):
		self.params = params
		
	# params is a dictionary of url parameters
	def setParams(self, params):
		self.params = params

	def setOneParam(self, param, value):
		self.params[param] = value

	def removeParam(self, param):
		if param in self.params:
			del self.params[param]

	def getImageURL(self):
		self.setOneParam("hgt.psOutput", "on")
		return self._getURL()
		
	def getLinkURL(self):
		self.removeParam("hgt.psOutput")
		return self._getURL()

	def _getURL(self):
		paramsenc = urllib.urlencode(self.params)
		url = "%s%s?%s" % (baseurl, pagename, paramsenc)
		return url

	def dumpPDF(self, pdffile):
		url = self.getImageURL()
		imgpagestring = urllib.urlopen(url).read()
		self._getPDFFromImagePage(imgpagestring, pdffile)

	def dumpPDFAndPNG(self, filebase):
		pdffile = filebase + ".pdf"
		self.dumpPDF(pdffile)
		pngfile = self.convertImage(pdffile, "png")
		return pdffile, pngfile

	# do we need both of these parameters?
	def setSessionDetails(self, sessionurl):
		self.setOneParam("hgS_doLoadUrl", "submit")
		self.setOneParam("hgS_loadUrlName", sessionurl)

	# in this case, bedurl could be a single bed file or
	# a url with a list of urls containing the bed files we want to upload
	def setBEDFiles(self, bedurl):
		self.setOneParam("hgt.customText", bedurl)

	def convertImage(self, inpath, outtype):
		# I wish python had extension removal in basename...
		basepath = inpath.replace(".pdf", "")
		outpath = basepath + "." + outtype
		# magic density value gives us the right appearance for a web output
		cmd = "%s -density 125 %s %s" % (convertcmd, inpath, outpath)
		print "executing convert command:", cmd
		os.system(cmd)
		return outpath

	# private
	def _getPDFLinkFromImagePage(self, pagestring):
		pdfre = 'HREF="([^"]+.pdf)"'
		match = re.search(pdfre, pagestring)
		return match.group(1)

	# private
	def _getPDFFromImagePage(self, pagestring, pdffile):
		link = self._getPDFLinkFromImagePage(pagestring)
		url = baseurl + link
		pdfdata = urllib.urlopen(url).read()
		pdffh = open(pdffile, "w")
		pdffh.write(pdfdata)

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: ./ucscimagegrab.py <sessionfile> <bedfileurl>"
		sys.exit(1)
	# these should both be publically available urls
	sessionurl = sys.argv[1]
	bedurl = sys.argv[2]
	params = { 	"hgt.psOutput"	:	"on",
				}
	ut = UCSCTools(params)
	ut.setSessionDetails(sessionurl)
	ut.setBEDFiles(bedurl)
	print ut.getLinkURL()
#	ut.dumpPDFAndPNG("tester")
	
	


