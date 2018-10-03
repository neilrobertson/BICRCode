#!/usr/bin/perl

use strict;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');

my $ga = $db->get_GeneAdaptor();

my $filename = "/home/pzs/codetesters/pytesters/exongeneids.txt";

open(FH, $filename);
my @genenames = ();

while(<FH>)
{
	chomp($_);
	if($_){ push(@genenames, $_); }
}

close(FH);

foreach my $genename (@genenames)
{
#	print "working at gene: ", $genename, "\n";
	my $gene = $ga->fetch_by_stable_id($genename);
	if($gene)
	{
		print $genename, ",", $gene->strand, ",", "chr", $gene->seq_region_name, ":", $gene->start, "-", $gene->end, "\n"
	}
	else
	{
		print STDERR "could not find: ", $genename, "\n";
	}
}

