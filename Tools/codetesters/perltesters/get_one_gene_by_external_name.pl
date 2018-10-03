#!/usr/bin/env perl

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);

my $ga = $registry->get_adaptor("human", "core", "gene");

my $geneid = $ARGV[0];

#print "looking for gene: " . $geneid . "\n";
my @genes = @ { $ga->fetch_all_by_external_name($geneid) };

#print "found: " . scalar(@genes) . " genes" . "\n";

foreach my $gene (@genes)
{
	print $geneid . "\t" . $gene->stable_id() . "\n";
}
