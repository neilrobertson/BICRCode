import unittest
import sys
import os

sys.path.append(os.path.expanduser("~/mount/repository/shared/baselib"))

from bed.bedToVstep import *

FILTER=''.join([(len(repr(chr(x)))==3) and chr(x) or '.' for x in range(256)])

def dump(src, length=8):
    N=0; result=''
    while src:
       s,src = src[:length],src[length:]
       hexa = ' '.join(["%02X"%ord(x) for x in s])
       s = s.translate(FILTER)
       result += "%04X   %-*s   %s\n" % (N, length*3, hexa, s)
       N+=length
    return result

class TestSequenceFunctions(unittest.TestCase):
    def compare(self, vstepWidth, withZeros, cores , input, expected):
        inputDataFilename = "unit-test-input-data.bed"
        inputData = open(inputDataFilename, "w")
        inputData.write(input)
        inputData.close()

        outputDataFilename = "unit-test-output-data.bed"
        doFile(inputDataFilename, outputDataFilename, vstepWidth, withZeros, cores)
        
        output = open(outputDataFilename, "U").read() # note using the U universal string format here in order to make the comparison work
        
        if expected.strip() == output.strip():
            success = True
        else:
            success = False
        
        os.remove(inputDataFilename)
        os.remove(outputDataFilename)
        
        if success == False:
            print
            print expected.strip()
            print "***"
            print output.strip()
            
        return success
    
    def testSimple(self):
        self.assertEqual(self.compare(200, False, 4,
'''chrY	1	150
chrY	201	300
chrY	350	400
chrY	401	600
chrY	650	750''', 
'''chrY	1	1
chrY	201	2
chrY	401	1
chrY	601	1'''), True)
        
	def testUnitBins(self):
		self.assertEqual(self.compare(1, False, 1,
'''chrY	1	20
chrY	5	15
chrY	12	23
chrY	22	25
chrY	27	28
''',
'''chrY	1	1
chrY	2	1
chrY	3	1
chrY	4	1
chrY	5	2
chrY	6	2
chrY	7	2
chrY	8	2
chrY	9	2
chrY	10	2
chrY	11	2
chrY	12	3
chrY	13	3
chrY	14	3
chrY	15	3
chrY	16	2
chrY	17	2
chrY	18	2
chrY	19	2
chrY	20	2
chrY	21	1
chrY	22	2
chrY	23	2
chrY	24	1
chrY	25	1
chrY	27	1
chrY	28	1'''), True)

    def testTinyBins(self):
        self.assertEqual(self.compare(10, False, 1,
'''chrY	1	150
chrY	205	300
chrY	350	400''',
'''chrY	1	1
chrY	11	1
chrY	21	1
chrY	31	1
chrY	41	1
chrY	51	1
chrY	61	1
chrY	71	1
chrY	81	1
chrY	91	1
chrY	101	1
chrY	111	1
chrY	121	1
chrY	131	1
chrY	141	1
chrY	201	1
chrY	211	1
chrY	221	1
chrY	231	1
chrY	241	1
chrY	251	1
chrY	261	1
chrY	271	1
chrY	281	1
chrY	291	1
chrY	341	1
chrY	351	1
chrY	361	1
chrY	371	1
chrY	381	1
chrY	391	1'''), True)

    def testHugeBins(self):
        self.assertEqual(self.compare(1000, False, 4,
'''chrY	1	150
chrY	201	300
chrY	350	400
chrY	401	600
chrY	650	750''', 
'''chrY	1	5'''), True)

    def testOverlapping(self):
        self.assertEqual(self.compare(200, False, 4,
'''chrY	1	150
chrY	101	200
chrY	125	160
chrY	140	400''',
'''chrY	1	4
chrY	201	1'''),  True)

    def testOverlappingTinyBins(self):
        self.assertEqual(self.compare(10, False, 1,
'''chrY	1	150
chrY	101	200
chrY	125	160
chrY	140	250''',
'''chrY	1	1
chrY	11	1
chrY	21	1
chrY	31	1
chrY	41	1
chrY	51	1
chrY	61	1
chrY	71	1
chrY	81	1
chrY	91	1
chrY	101	2
chrY	111	2
chrY	121	3
chrY	131	4
chrY	141	4
chrY	151	3
chrY	161	2
chrY	171	2
chrY	181	2
chrY	191	2
chrY	201	1
chrY	211	1
chrY	221	1
chrY	231	1
chrY	241	1'''),True)

    def testTouchingIntervals(self):
        self.assertEqual(self.compare(200, False, 4,
'''chrY	100	200
chrY	200	210
chrY	210	310
chrY	310	330''', 
'''chrY	1	2
chrY	201	3'''), True)

    def testDuplicateIntervals(self):
        self.assertEqual(self.compare(200, False, 4,
'''chrY	100	200
chrY	100	200
chrY	100	200
chrY	210	310
chrY	210	310
chrY	210	310
chrY	210	310''', 
'''chrY	1	3
chrY	201	4'''), True)

    def testBigGaps(self):
        self.assertEqual(self.compare(200, False, 4,
'''chrY	100	200
chrY	150	250
chrY	1500	1610
chrY	11000	11090''', 
'''chrY	1	2
chrY	201	1
chrY	1401	1
chrY	1601	1
chrY	10801	1
chrY	11001	1'''), True)

    def testUnsortedInput(self):
        self.assertEqual(self.compare(200, False, 4,
'''chrY	401	600
chrY	350	400
chrY	1	150
chrY	201	300
chrY	650	750''', 
'''chrY	1	1
chrY	201	2
chrY	401	1
chrY	601	1'''), True)

    def testIntervalsWithinIntervals(self):
        self.assertEqual(self.compare(200, False, 4,
'''chrY	100	1000
chrY	200	250
chrY	225	275
chrY	500	600
chrY	525	625
chrY	650	700''', 
'''chrY	1	2
chrY	201	3
chrY	401	3
chrY	601	3
chrY	801	1'''), True)

if __name__ == "__main__":
    unittest.main()
