#!/usr/bin/env perl

use Bio::EnsEMBL::Registry;

## connect to database
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);


my $species = "human";
my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );

my @chromosomes = (	
					"NT_113953"
#					"NT_113944"
#					"NT_113933"
#					"NT_113925"
#					"NT_113886"
#					"NT_113871"
# 					"1", 
#  					"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
#  					"12", "13", "14", "15", "16", "17", "18", "19",  "X", "Y" 
					);

#print $slice_adaptor->dbc()->dbname . "\n";

foreach $chr (@chromosomes) {           # make an array of chromosome names
	my $chr_slice = $slice_adaptor->fetch_by_region($coord_system_name,$chr) or die "Failed to slice chromosome $chr";
	
	print join("\t", ( $chr, $chr_slice->end )) . "\n";
}
