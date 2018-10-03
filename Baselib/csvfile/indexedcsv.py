import csv
import sys
import collections
import ordereddict # sudo pip --proxy http://wwwcache.gla.ac.uk:8080 install ordereddict


class IndexedCSV(ordereddict.OrderedDict):
    def __init__(self, csvfileloc, header=True, keyPos = 0, key=None, useRowNumAsKey=False, defaultkeys=[]):
        
        super(IndexedCSV, self).__init__()
        #self.clear()
        
        self.keys = defaultkeys
        valuesInputFile = csv.reader(open(csvfileloc, "r"), delimiter='\t')
        rowsRead = 0
        
        self.rows = {}
        
        for row in valuesInputFile:
            if len(row)==0:
                continue
            elif row[0].startswith("#"):
                self.keys = []
                for k in row:
                    self.keys.append(k.replace("#",""))
                #print >> sys.stderr, "Loading Indexed CSV file with keys: "+repr(self.keys)
                header = False
                if key!=None:
                    keyPos = self.keys.index(key)
                
            elif header==True:
                # first column is keys/headers
                self.keys = []
                for k in row:
                    self.keys.append(k.replace("#",""))
                header=False # dealt with header
                if key!=None:
                    keyPos = self.keys.index(key)
                #print >> sys.stderr, "Loading Indexed CSV file with keys: "+repr(self.keys)
            else:
                rowsRead+=1
                # actual data                
                
                data = {}
                for i in range(0, len(self.keys)):
                    data[self.keys[i]]=row[i]                
                
                keyValue = rowsRead if useRowNumAsKey else row[keyPos]
                if keyValue in self:
                    assert "Seen "+keyValue+" already!"
                
                self.rows[keyValue] = row
                self[keyValue] = data

        assert rowsRead == len(self), "Number of rows read != number of entries in indexed csv -- Possible Duplicate keys? " + str(rowsRead) + "/" + str(len(self))

class ColumnIndex(collections.defaultdict):
    def __init__(self,indexedcsv,column):
        
        # call the super constructor for default dict and make it a list
        super(ColumnIndex,self).__init__(list)
        
        for row in indexedcsv:
            if column in indexedcsv[row]:
                self[indexedcsv[row][column]].append(row)
            else:
                self[""].append(row)