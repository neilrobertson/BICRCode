"""
Utilities for dealing with the TFD motif collection.
"""

import os.path
from motility.objects import PWM

class TFD(object):

    def load(self,list_file, matrix_dir):
        """
        Load all of the matrices in list_file from matrix_dir.  Return a
        dictionary of PWMs.
        """
    
        l = open(list_file).readlines()
        
        l = [ i.strip() for i in l ]
    
        d = {}
        for filename in list:
            name = filename.split('.')[0]
            filename = os.path.join(matrix_dir, filename)
    
            lines = open(filename).readlines()
            
            lines2 = []
            for l in lines:
                if l[0:2] == 'MA' and l[5:8] == '/M:':
                    assert l[9:13] == 'M=  '
                    l = l[13:].split(';')[0]
    
                    A, C, G, T = map(int, l.split(','))
                    lines2.append((A, C, G, T, 0))
    
            d[name] = PWM(lines2)
    
        return d
