#!/usr/bin/env perl

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);

my $ga = $registry->get_adaptor("human", "core", "gene");

my $gene = $ga->fetch_by_stable_id("ENSG00000012048");

print $gene->stable_id() . " " . $gene->start() . " " . $gene->end() . "\n";
