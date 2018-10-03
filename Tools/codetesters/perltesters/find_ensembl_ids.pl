#!/usr/bin/env perl

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $filename = "/home/pzs/workspace/biomart/ensemblids.txt";

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');
my $ga = $db->get_GeneAdaptor();

open(IN, $filename);

while (<IN>) {
	$genename = $_;
	chomp($genename);
	
	$gene = $ga->fetch_by_stable_id($genename);
	if($gene)
	{
		print "found gene name ", $genename, "\n";
		print $gene->strand, "\n";
	}
	else
	{
		print "could not find gene ", $genename, "\n";
	}
	
}

close(IN);
