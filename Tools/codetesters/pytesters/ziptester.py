#!/usr/bin/env python

import zipfile
import os

fname = "pythonzip.zip"

directory = "svntester"
zfile = zipfile.ZipFile(fname, "w")
def walker( zip, directory, files, root=directory ):
	for file in files:
		file = os.path.join( directory, file )
		# yes, the +1 is hacky...
		archiveName = file[len(os.path.commonprefix( (root, file) ))+1:]
		zip.write( file, archiveName, zipfile.ZIP_DEFLATED )
		print file
os.path.walk( directory, walker, zfile  )
zfile.close()


