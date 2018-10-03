'''
Created on 5 Jan 2015

@author: neilrobertson
'''



class FileConfiguration(object):
    def __init__(self, *args, **kwargs):
        pass
    
    
    def test_HiC_bedFile(self, inputFilename):
        
        with open(inputFilename, "r") as inputFile:
            pass