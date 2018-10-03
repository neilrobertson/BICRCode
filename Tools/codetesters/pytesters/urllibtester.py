#!/usr/bin/env python

import urllib

blast_wsdUrl='http://www.ebi.ac.uk/Tools/webservices/wsdl/WSNCBIBlast.wsdl'
fh = urllib.urlopen(blast_wsdUrl)

wsdlfile = fh.read()

print wsdlfile

