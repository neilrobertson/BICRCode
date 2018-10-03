class Column(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value
    
    def getName(self):
        return self.name
    
    def toString(self):
        return str(self.value)

class ListColumn(Column):
    def __init__(self, name, values):
        super(ListColumn, self).__init__(name, values)
    
    def toString(self):
        return ", ".join(self.value)

class Classifiee():
    def __init__(self,  chr, start, stop):
        self.chr = chr
        self.start = start
        self.stop = stop
        assert self.start <= self.stop,  "Stop > start"
        self.columns = dict()
        
        self.addColumn(Column("chr", chr))
        self.addColumn(Column("start", start))
        self.addColumn(Column("stop", stop))
    
    def addColumn(self, column):
        self.columns[column.getName()] = column
        
    def getRawColumnValue(self,columnname):
        return self.columns[columnname].value
    
    def getColumnValue(self,columnname):
        column = self.columns[columnname]
        if column.value == True:
            return 1
        elif column.value == False:
            return 0
        else:
            return str(column.toString())
