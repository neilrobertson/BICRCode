#!/usr/bin/perl
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);  
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db
(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

$slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "slice" );
# $my $slice_adaptor = $registry->get_adaptor('Human','Core','Slice');


@chromosomes = (	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
					"12", "13", "14", "15", "16", "17", "18", "19", "20", "21", 
					"22", "X", "Y" );

my $numgenes = 0;

foreach my $chr (@chromosomes) {           # make an array of chromosome names
	print "working at chromosome: ". $chr . "... ";

#Inside the chromosome loop, use the slice adaptor to get the whole chromosome, like this:
my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr);

#inside the chromosome loop, use the gene adaptor to get all the genes from this slice into a Perl list:

my @genes = @{$slice->get_all_Genes()};

$numgenes += scalar(@genes);

print scalar(@genes) . "\n";

# print "id","\t"; 
# print "region","\t"; 
# print "start","\t"; 
# print "end","\n"; 
# 
# foreach  $genes(@genes) {
# 	print $genes->display_id(), "\t";
# 	print $genes->seq_region_name(), "\t";
# 	print $genes->start(), "\t";
# 	print $genes->end(), "\n";
# }
# 

}

print $numgenes . "\n";

# Loop over the gene objects in the list. For each gene, print the display_id, the chromosome (use seq_region_name), the start and the end for each gene.

