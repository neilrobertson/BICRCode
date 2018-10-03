#!/usr/bin/env perl

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $filename = "/home/pzs/workspace/biomart/vegaids.txt";

my $vdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(	-host => 'ensembldb.ensembl.org',
												-user => 'anonymous',
												-dbname => 'homo_sapiens_vega_36_35i');

my $va = $vdb->get_GeneAdaptor();

open(IN, $filename);

while (<IN>) {
	$genename = $_;
	chomp($genename);
	
	$gene = $va->fetch_by_stable_id($genename);
	if($gene)
	{
		print "found gene name ", $genename, "\n";
	}
	else
	{
		print "could not find gene ", $genename, "\n";
	}
	
}

close(IN);
