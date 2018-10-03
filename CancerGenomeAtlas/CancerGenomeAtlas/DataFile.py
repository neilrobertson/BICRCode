'''
Created on 1 Sep 2014

@author: neilrobertson
'''

class DataFile(object):
    def __init__(self, fileName, columnDelimiter = "\t"): 
        self.fileName = None
              
        
    def openFileAndParseHeaders(self):
        ''' 
        Opens the RNA data file in TCGA format and returns header values
        '''
        assert False, "Open file and parse headers method not defined"
    
    def parseValues(self):
        '''
        Parses all data lines in file and yields them as a generator
        '''
        assert False, "Parsing not defined"
                
    def closeFile(self):
        '''
        Closes data file
        '''
        assert False, "Close file not defiined"