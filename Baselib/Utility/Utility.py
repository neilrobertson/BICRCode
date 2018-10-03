'''
Created on 2 Sep 2014

A series of static methods (or otherwise) that can help with routine pythonic problems.

@author: neilrobertson
'''

class Utility(object):
    def __init__(self, *args, **kwargs):
        pass
    
    @staticmethod
    def getColumns(inFile, delim="\t", header=True):
        """
        Get columns of data from inFile. The order of the rows is respected
        
        :param inFile: column file separated by delim
        :param header: if True the first line will be considered a header line
        :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that 
        are headings in the inFile, and values are a list of all the entries in that
        column. indexToName dict maps column index to names that are used as keys in 
        the cols dict. The names are the same as the headings used in inFile. If
        header is False, then column indices (starting from 0) are used for the 
        heading names (i.e. the keys in the cols dict)
        """
        cols = {}
        indexToName = {}
        for lineNum, line in enumerate(inFile):
            if lineNum == 0:
                headings = line.split(delim)
                i = 0
                for heading in headings:
                    heading = heading.strip()
                    if header:
                        cols[heading] = []
                        indexToName[i] = heading
                    else:
                        # in this case the heading is actually just a cell
                        cols[i] = [heading]
                        indexToName[i] = i
                    i += 1
            else:
                cells = line.split(delim)
                i = 0
                for cell in cells:
                    cell = cell.strip()
                    cols[indexToName[i]] += [cell]
                    i += 1 
        return cols, indexToName
    
    @staticmethod
    def unpack_tuple(inputTupleOrList):
        for item in inputTupleOrList:
            yield item
            
    @staticmethod
    def unpack_fromPositions(seq, positions):
        for pos in positions:
            yield seq[pos]

    @staticmethod
    def unpack_nfirst(seq, nfirst):
        it = iter(seq)
        for x in xrange(nfirst):  # @UnusedVariable
            yield next(it, None)
        yield tuple(it)
        
    @staticmethod
    def checkDirectory(directory):
        import os
        if os.path.isdir(directory):
            if directory[-1] != r"/":
                directory = directory + r"/"
            return directory
        else:
            return False
        
    @staticmethod
    def create_directory(directory):
        import os
        if not os.path.exists(directory): 
            try:
                os.makedirs(directory)
                return True
            except:
                return False
            return True