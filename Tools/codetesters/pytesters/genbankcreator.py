#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Alphabet import IUPAC
import sys

from bp import breakpoint

sf = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(3,10), type="misc_feature", qualifiers = { "label" : [ "testlabel" ]})

sr = SeqRecord(Seq("GTTACGACATCGACTAGAGGCATAGCA", IUPAC.unambiguous_dna), 
				id="tester", name="testertester", description="first stab",
				features=[sf])

SeqIO.write([sr], sys.stdout, "genbank")
