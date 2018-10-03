#!/usr/bin/env perl

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $filename = "/home/pzs/workspace/biomart/plainnames.txt";

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');
my $ga = $db->get_GeneAdaptor();

open(IN, $filename);

while (<IN>) {
	$genename = $_;
	chomp($genename);
	
	@genes = @{$ga->fetch_all_by_external_name($genename)};
	if(scalar(@genes > 0))
	{
		print "found gene name ", $genename, "\n";
	}
	else
	{
		print "could not find gene ", $genename, "\n";
	}
	
}

close(IN);
