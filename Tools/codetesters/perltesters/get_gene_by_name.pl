#!/usr/bin/env perl

use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);

my $ga = $registry->get_adaptor("human", "core", "gene");

my $genename = $ARGV[0];

my @genes = @{$ga->fetch_all_by_external_name($genename)};
if(scalar(@genes) == 1)
{
	my $gene = $genes[0];
	print $genename, ",", $gene->strand, ",", "chr", $gene->seq_region_name, ":", $gene->start, "-", $gene->end, "\n"
}

