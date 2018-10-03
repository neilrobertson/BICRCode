#!/usr/bin/env python

class ExecTester:
	def __init__(self):
		exec("self.testvar = 3")
		
	def testfunc(self):
		print self.testvar
		
et = ExecTester()
et.testfunc()
