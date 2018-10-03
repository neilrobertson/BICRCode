import csv
import os
import sys
# sudo pip install pyExcelerator
from pyExcelerator import *

class CSVtoExcelWriter():
    def __init__(self, csvreader, outputfileLoc):
        
        wb = Workbook()
        sheet = wb.add_sheet('CSV File')
        
        x = 0
        for row in csvreader:
            for y in range(len(row)):
                # try to write it as a float first, else string
                try:
                    sheet.write(x,y,float(row[y]))
                except ValueError:
                    sheet.write(x,y,row[y])
                sheet.write(x,y,row[y])
            x += 1
        wb.save(outputfileLoc)
    
class ExcelFile():
    def __init__(self, outputfileLoc):
        self.outputfileLoc = outputfileLoc
        
        self.wb = Workbook()
        self.sheet = self.wb.add_sheet('CSV File')
        
        self.x = 0
    
    def writerow(self, row):
        for y in range(len(row)):
            self.sheet.write(self.x,y,row[y])
        self.x += 1
    
    def save(self):
        self.wb.save(self.outputfileLoc)
