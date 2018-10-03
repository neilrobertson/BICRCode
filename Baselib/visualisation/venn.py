#!/usr/bin/python
"""PyVenn

http://code.google.com/p/pyvenn/

Venn Diagrams for python,
both proportional and 'normal'.

Passes diagrams (normal) on to the R library 'vennerable'
using rpy2 (and numpy).


Usage:
    vd = VennDiagram ( [ ("One set", [1,2,3]), ("Other Set", [3, 4, 5, 6])])
    vd.plot_normal('test_normal.png')

Requires:
    numpy

Option:
    rpy2 ( and the Vennerable R library)

Released under:
The MIT License

Copyright (c) 2010, Florian Finkernagel

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


"""

import math
import numpy
import unittest
import rpy2.robjects
import rpy2.robjects.numpy2ri #so the numpy->r interface works
import rpy2.rinterface
import csv

rpy2.robjects.numpy2ri.activate()

_r_loaded = False
robjects = None
def load_r():
    global _r_loaded
    global robjects
    robjects = rpy2.robjects
    if not _r_loaded:
        robjects.r('library(Vennerable)')

class VennDiagram:

    def __init__(self, name_set_tuples_or_dict):
        if hasattr(name_set_tuples_or_dict, 'items'):
            sets = name_set_tuples_or_dict.items()
        else:
            sets = name_set_tuples_or_dict
        self.sets = []
        for name, group in sets:
            self.sets.append((name, set(group)))

    default_colors = [ (228 / 255.0, 26 / 255.0, 28 / 255.0),
                              (55 / 255.0, 126 / 255.0, 184 / 255.0),
                              (77 / 255.0, 175 / 255.0, 74 / 255.0)]

    def plot_normal(self, output_filename, width=8):
        self._venn_plot_weights(output_filename, width, width)

    def getFriendlyName(self,string):
        return string[string.rfind("/")+1:]

    def _get_set_names(self):
        return [self.getFriendlyName(x[0]) for x in self.sets]

    def _get_set_values(self):
        return [x[1] for x in self.sets]

    def _venn_plot_sets(self, output_filename, width=8, height=8):
        """Plot a venn diagram into the pdf file output_filename.
        Takes a dictionary of sets and passes them straight on to R"""
        #raise TypeError("This function should no longer be used. use _venn_plot_weigths instead")
        sets = self.sets
        load_r()
        
        if output_filename.endswith('.pdf'):
            type = "pdf"
        elif output_filename.endswith('.png'):
            type = "png"
            dpi = 72 # pixels per inch
            width = width * dpi
            height = height * dpi
        else:
            raise ValueError("_venn_plot_weights currently only understands .pdf and .png output filenames")
        
        robjects.r(type)(output_filename, width=width, height=height)
        
        x = robjects.r('Venn')(Sets = [numpy.array(list(x)) for x in self._get_set_values()], SetNames=[x for x in self._get_set_names()])
        
        if len(self.sets) <= 3:
            venn_type = 'circles'
            weights = True
            euler = True
        elif len(self.sets) == 4:
            venn_type = 'squares'
            weights = False
            euler = False
        else:
            venn_type = 'ChowRuskey';
            weights = False
            euler = False
        
        robjects.r('plot')(x, **{'type': venn_type, 'doWeights': weights, 'doEuler':euler})
        robjects.r('dev.off()')


    def _venn_plot_weights(self, output_filename, width=8, height=8):
        """Plot a venn diagram into the pdf file output_filename.
        Takes a dictionary of sets and does the intersection calculation in python
        (which hopefully is a lot faster than passing 10k set elements to R)
        """
        load_r()
        weights = [0]
        sets_by_power_of_two = {}
        for ii, kset in enumerate(self._get_set_names()):
            iset = self._get_set_values()[ii]
            sets_by_power_of_two[2**ii] = set(iset)
        for i in xrange(1, 2**len(self.sets)):
            sets_to_intersect = []
            to_exclude = set()
            for ii in xrange(0, len(self.sets)):
                if (i & (2**ii)):
                    sets_to_intersect.append(sets_by_power_of_two[i & (2**ii)])
                else:
                    to_exclude = to_exclude.union(sets_by_power_of_two[(2**ii)])
            final = set.intersection(*sets_to_intersect) - to_exclude
            weights.append( len(final))
        
        if output_filename.endswith('.pdf'):
            type = "pdf"
        elif output_filename.endswith('.png'):
            type = "png"
            dpi = 72 # pixels per inch
            width = width * dpi
            height = height * dpi
        else:
            raise ValueError("_venn_plot_weights currently only understands .pdf and .png output filenames")
        
        robjects.r(type)(output_filename, width=width, height=height)
        x = robjects.r('Venn')(Weight = numpy.array(weights), SetNames=self._get_set_names())
        
        if len(self.sets) <= 3:
            venn_type = 'circles'
            weights = True
            euler = True
        elif len(self.sets) == 4:
            venn_type = 'squares'
            weights = False
            euler = False
        else:
            venn_type = 'AWFE'
            weights = False
            euler = False

        try:
            print venn_type
            print x
            robjects.r('plot')(x, **{'type': venn_type, 'doWeights': weights, 'doEuler': euler})
        except rpy2.rinterface.RRuntimeError:
            # sometimes the weighted versions will fail if it is undrawable
            robjects.r('plot')(x, **{'type': venn_type, 'doWeights': False, 'doEuler': False})
         
        robjects.r('dev.off()')

class TestVennDiagram(unittest.TestCase):
    """These 'unit tests' don't deserve their name - but at least the
    ascertain that we don't throw an exception."""
    def test_normal(self):
        of = "test.png"
        sets = {
            'A': [1,2,3,44],
            'B': [2,3,55,66,77]
        }
        d = VennDiagram(sets)
        d.plot_normal(of)

    

    

if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    def printUsage():
        print "venn.py -f output.png (--width 8) a b c d etc"
        sys.exit(1)
    parser = OptionParser()
    parser.add_option('-t', '--test', dest="is_test", action="store_true", help="Run unittests (and create a bunch of test*.png/pdf")
    parser.add_option('-f', '--output_filename', dest='output_filename', help="output filename, must end with .png or .pdf")
    parser.add_option('-w', '--width', dest='width', help="output image width in inch (resolution: 72 dpi")
    
    (options, args) = parser.parse_args()
    if options.is_test is True:
        unittest.main()
    if options.width is None:
        width = 8
    else:
        width = int(options.width)
    if options.output_filename is None:
        printUsage()

    sets = list()
    
    for arg in args:
        csvreader = csv.reader(open(arg,"r"))
        s = set()
        for row in csvreader:
            s.add(row[0].strip())
        sets.append((arg,s))
    
    vd = VennDiagram(sets)
    vd.plot_normal(options.output_filename, width)
    sys.exit(0)
