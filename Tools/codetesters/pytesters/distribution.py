#!/usr/bin/env python

from pylab import randn, hist, show
from numpy.random import *

import pdb
_pdb = pdb.Pdb()
breakpoint = _pdb.set_trace


x = uniform(size=100000)
breakpoint()
hist(x, 100)
show()
