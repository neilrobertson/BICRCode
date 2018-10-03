#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
sequence = Seq("AAA", generic_dna)
sequencerecord = SeqRecord(sequence, id="testid", description="This is just a test")
output_handle = open("example.fasta", "w")
SeqIO.write([ sequencerecord ], output_handle, "fasta")
output_handle.close()
