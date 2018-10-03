import sys
import os
from datastructures.binarytree import *

class GenomeBinaryTree(object):
    def __init__(self):
        self.values = {}
    
    def insert(self, chr, key, value):
        if not chr in self.values:
            self.values[chr] = binary_tree()
        self.values[chr][key] =value
    
    def hasChromosome(self, chr):
        return chr in self.values
    
    def getChromosomeBinaryTree(self, chr):
        if self.hasChromosome(chr):
            return self.values[chr]
    
    def getNearestNodes(self, chr, key):
        # get the interval tree for here
        if not self.hasChromosome(chr):
            return [] # its on a chromosone we dont have any values for
        tree = self.getChromosomeBinaryTree(chr)
        
        return tree.findnearest(key)
    
    def getNearestValues(self, chr, key):
        nodes = self.getNearestNodes(chr, key)
        
        values = []
        for n in nodes:
            values.append(n.data)
            
        return values
