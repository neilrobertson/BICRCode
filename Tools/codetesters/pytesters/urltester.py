#!/usr/bin/env python

import urllib

f = urllib.urlopen("http://scholar.google.co.uk/scholar?q=schoeberl&hl=en&lr=&btnG=Search")
print f.read()
