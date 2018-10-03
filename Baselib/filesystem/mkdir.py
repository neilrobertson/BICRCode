'''
Created on 23 Aug 2010

@author: mcbryan
'''
import os

def makeDirectory(dir):
    if not os.access(dir, os.R_OK):
        os.mkdir(dir)
        return dir