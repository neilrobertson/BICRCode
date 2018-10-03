#!/usr/bin/env python

# Import WSDL package
from SOAPpy import WSDL
import os
import sys
 
# Create service interface
wsdlUrl = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSFasta.wsdl'
server = WSDL.Proxy(wsdlUrl)
 
# Create the service interface
# Configure HTTP proxy from OS environment (e.g. http_proxy="http://proxy.example.com:8080")
if os.environ.has_key('http_proxy'):
	http_proxy_conf = os.environ['http_proxy'].replace('http://', '')
elif os.environ.has_key('HTTP_PROXY'):
	http_proxy_conf = os.environ['HTTP_PROXY'].replace('http://', '')
else:
	http_proxy_conf = None
server.soapproxy.http_proxy = http_proxy_conf

if len(sys.argv) != 2:
	print "wrong!"
	sys.exit(1)
jobid = sys.argv[1]

print dir(server.getResults(jobid)[0])
