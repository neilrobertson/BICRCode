#!/usr/bin/env python

from ensembl import *

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace

serverRegistry = get_registry(host='ensembldb.ensembl.org', user='anonymous')
