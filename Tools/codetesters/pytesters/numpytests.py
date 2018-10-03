#!/usr/bin/env python

from numpy import *

def nangenerator(i, j):
	return i + j + nan

def nangenerator1(i):
	return i + nan

print fromfunction(nangenerator, (1, 10))[0]
