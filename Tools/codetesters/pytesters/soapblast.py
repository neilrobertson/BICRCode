#!/usr/bin/env python

# Import WSDL package
from SOAPpy import WSDL
import os
 
# Create service interface
wsdlUrl = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSWUBlast.wsdl'
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

 
# Query sequence
seq = """>Q8E5Q5_STRA3
MKLSKRYRFWQKVIKALGVLALIATLVLVVYLYKLGILNDSNELKDLVHKYEFWGPMIFI
VAQIVQIVFPVIPGGVTTVAGFLIFGPTLGFIYNYIGIIIGSVILFWLVKFYGRKFVLLF
MDQKTFDKYESKLETSGYEKFFIFCMASPISPADIMVMITGLSNMSIKRFVTIIMITKPI
SIIGYSYLWIYGGDILKNFLN"""
 
# Structure containing parameters
blast_params = {
  'program':'blastp',      # Protein vs. protein search
  'database':'swissprot',  # Database to search
  'email':'pzs@dcs.gla.ac.uk',    # User e-mail address
  'async':1                # Async submission
}
 
# Structure containing the sequence data
blast_data = [{'type':'sequence', 'content':seq}]
 
# Submit the job to the service passing the data structures
jobId = server.runWUBlast(params=blast_params,content=blast_data)
print jobId
