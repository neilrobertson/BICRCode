'''
Created on 21 Sep 2010

@author: mcbryan
'''
import os
import csv

class ChrList(list):
    def __init__(self,build):
        infile = os.path.expanduser("~/mount/publicdata/"+build+"/chrList."+build)
        
        reader = csv.reader(open(infile), delimiter="\t")
        
        for row in reader:
            self.append(row[0])