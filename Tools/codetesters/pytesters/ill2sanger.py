#!/usr/bin/env python

from Bio import SeqIO

from Bio import SeqIO
SeqIO.write(SeqIO.parse(open("/home/pzs/solexa/run1/s_4_sequence.txt", "rU"), "fastq-illumina"), \
            open("tester.fq", "w"), "fastq")
