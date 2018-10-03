#!/usr/bin/perl

use strict;
use warnings;
use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");
use Bio::SeqIO;
use Bio::EnsEMBL::Registry;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');

my $sa = $db->get_SliceAdaptor();



my @chromosomes = (	
					"1", 
 					"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
 					"12", "13", "14", "15", "16", "17", "18", "19",  
 					"20", "21", "22",
					"X", "Y" 
					);

foreach my $chr (@chromosomes) {           # make an array of chromosome names
	my $slice = $sa->fetch_by_region('chromosome',$chr) or die "Failed to slice chromosome $chr";
	
	my @genes = @{ $slice->get_all_Genes() };
	foreach my $gene (@genes) {
		my $geneid = $gene->stable_id();
		my $genename = $gene->external_name();
		my $chrom = $gene->seq_region_name();
 		my $start = $gene->start();
 		my $end = $gene->end();
		if (!$genename) { $genename = ""; }
		print join("\t", ($geneid, $genename, $chrom, $start, $end)) . "\n"
	}	
}



