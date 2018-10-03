#!/usr/bin/perl -w


use strict;
use warnings;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;




## connect to database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');


## get a slice adaptor
my $slice_adaptor = $db->get_SliceAdaptor();

my $chromcount = 0;

## foreach chromosome in turn
foreach my $chr (1..22, "X", "Y") {
	print "working at chromosome ", $chr, "\n";

	## get the chromosome object
	my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr);
	
	my @genes = @{$slice->get_all_Genes_by_type('protein_coding', 'ensembl')};
	$chromcount += scalar(@genes);
	## for all genes on the chromosome
# 	foreach my $gene (@{$slice->get_all_Genes_by_type('protein_coding', 'ensembl')}) {
# 
# 		my $genename = $gene->display_id();
# 		my $geneextname = $gene->external_name();
# 		if(!$geneextname) { $geneextname = ""; }
# 		my $chrom = "chr" . $gene->seq_region_name();
# 		my $start = $gene->start();
# 		my $end = $gene->end();
# 		my $strand = $gene->strand();
# 		my @line = ($geneextname, $genename, $chrom, $start, $end, $strand);
# 		print  join("\t", @line, "\n");
# 
# 	}
}

print "found ", $chromcount, "\n";
