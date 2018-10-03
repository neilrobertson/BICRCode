#!/usr/bin/perl

use strict;
use warnings;
use lib ("/home/pzs/src/ensembl-54/ensembl/modules/");
use Bio::SeqIO;

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";

$reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');

my $species = "human";

my $sa = $reg->get_adaptor($species, 'core', 'Slice');



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
		print join("\t", ($gene->stable_id(), $gene->external_name(), $gene->strand(), $gene->seq_region_name(), $gene->start(), $gene->end())) . "\n"
	}	
}



