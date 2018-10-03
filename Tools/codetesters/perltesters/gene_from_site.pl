#!/usr/bin/env perl

use Bio::EnsEMBL::Registry;

## connect to database
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);


my $species = "human";
if(scalar(@ARGV) != 3) {
	die "wrong number of arguments!";
}

my $chr = $ARGV[0];
my $start = $ARGV[1];
my $end =   $ARGV[2];

my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );

my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $start, $end) or die "Failed to slice chromosome $chr";

@genes = @{ $slice->get_all_Genes() };

print "found: " . scalar(@genes) . " genes\n";

foreach $gene ( @genes ) {
	print $gene->stable_id() . "\n";
	print $gene->external_name() . "\n";
}

