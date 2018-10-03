'''
Created on 8 Jan 2014

@author: mcbryan
'''

import getopt
import sys

import urllib2
import re

# sudo pip install poster
from poster.encode import multipart_encode
from poster.streaminghttp import register_openers

register_openers()

def UTRelements(fasta):
    regRNAbase = "http://regrna2.mbc.nctu.edu.tw/"
    regRNAexec = "detection_output.php"
    
    values = {'S1' : fasta,
              'UTRsite' : 'ON' }
    
    data, headers = multipart_encode(values)
    
    headers.update({ 'User-Agent' : 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:26.0) Gecko/20100101 Firefox/26.0',
                'Accept' : 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8', 
                'Accept-Language' : 'en-gb,en-us;q=0.5'})
    
    #data = urllib.urlencode(values)
    req = urllib2.Request(regRNAbase + regRNAexec, data, headers)
    response = urllib2.urlopen(req)
    
    results = None
    for line in response.readlines():

        if "Error" in line:
            raise Exception("Error in input data")

        if "Tab-Delimited" in line:
            results = re.sub('.*Results','Results',re.sub('. target=._.*', '', line))
    
    UTRelements = set()
    
    if results != None:
        req = urllib2.Request(regRNAbase + results)
        response = urllib2.urlopen(req)
        
        for line in response.readlines()[1:]:
            UTRelements.add(line.split("\t")[1].strip())
    
    return UTRelements


if __name__ == '__main__':

    fasta = '''>X83878|B.subtilis xpt and pbuX genes|operon; xanthine permease; xanthine phosphoribosyltransferase
aatattcaaatctctatctgttataatcaaaagcctggcggcgcggtcgtcagactcttt
tatatcgaatccccttgaaatacgaatgatatctaaaaaaacaaaattaaagttcgggaa
tttttattttcagcctatgcaagagattagaatcttgatataatttattacaatataata
ggaacactcatataatcgcgtggatatggcacgcaagtttctaccgggcaccgtaaatgt
ccgactatgggtgagcaatggaaccgcacgtgtacggttttttgtgatatcagcattgct
tgctctttatttgagcgggcaatgctttttttattctcataacggaggtagacaggatgg
aagcactgaaacggaaaatagaggaagaaggcgtcgtattatctgatcaggtattgaaag
tggattcttttttgaatcaccaaattgatccgctgcttatgcagagaattggtgatgaat
ttgcgtctaggtttgcaaaagacggtattaccaaaattgtgacaatcgaatcatcaggta
tcgctcccgctgtaatgacgggcttgaagctgggtgtgccagttgtcttcgcgagaaagc
ataaatcgttaacactcaccgacaacttgctgacagcgtctgtttattcctttacgaagc
aaacagaaagccaaatcgcagtgtctgggacccacctgtcggatcaggatcatgtgctga
ttatcgatgattttttggcaaatggacaggcagcgcacgggcttgtgtcgattgtgaagc
aagcgggagcttctattgcgggaatcggcattgttattgaaaagtcatttcagccgggaa
gagatgaacttgtaaaactgggctaccgagtggaatctttggcaagaattcagtctttag
aagaaggaaaagtgtccttcgtacaggaggttcattcatgagaaatggattcggcaaaac
gctgtctttagggattcagcatgttcttgccatgtatgccggcgccattgtcgttcctct
gattgtcggaaaagcaatgggactgactgtcgagcagctgacttacttagtatcgattga
tatttttatgtgcggtgtggctacacttctgcaagtgtggagcaaccgattttttgggat
cgggcttccggtagtgcttggctgtacctttacagctgtatcgccgatgatagcgattgg
atctgaatatggggtttcaacagtttacggcagcattatcgcttcaggcattcttgtcat
tcttatttcatttttctttggaaagctcgtatcgttttttccgccggtcgtgacaggctc
tgttgttacgattatcggtattacactgatgccggttgccatgaataacatggccggcgg
agaaggaagtgcagatttcggagatctctccaatcttgcacttgcttttaccgtgttgag
tatcattgtgcttctataccgttttacaaaaggctttatcaagtccgtctcgattttgat
cggtattttgattggcaccttcatcgcatattttatgggaaaagttcaatttgataatgt
ttcggacgcggcagttgttcaaatgattcagccattttacttcggagcgccgtcttttca
cgcagcgcctatcattacgatgtccatcgttgcaattgtcagccttgtggagtcaactgg
tgtttactttgctttaggtgacctgacaaaccgccgtttgacagagatagatttgtctaa
aggctatcgagctgaaggactggctgtgcttctgggcggtatttttaacgcttttcctta
cacggcattctctcaaaatgtcgggcttgttcagctcacagggatcaaaaaaaatgccgt
cattgtggtcacaggtgtgatcctaatggcatttgggctgtttccaaagatcgctgcttt
cacgactattattccatccgctgttttgggcggtgcaatggtggcgatgtttggaatggt
gattgcttacggaattaaaatgctcagccgcattgattttgcgaaacaggaaaatctgct
gattgttgcttgctcagtgggattggggctcggtgtgaccgttgtgccggatattttcaa
acagcttccgagtgccttaacgctgcttacaacgaatggcattgtggcgggcagctttac
tgcagtcgtcttaaatattgtatataacgttttttctaaagctaaaaaaatagaacaaga
agctgacctcgctgaacaaaaaacagcagtctaactccgccgcggcggagttttttttgc
atataaaaaagattgttgcggcacaatgtatgcaaagaggtgatcgcatggcgtttattt
tatccattggaacaagtcttcctgcgtataatgtaaaccaagagaaagcggccgagttcg
cccgctatatgtt'''

    print UTRelements(fasta)

