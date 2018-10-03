#!/usr/bin/perl

use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);

my $sa = $registry->get_adaptor("human", "core", "slice");

my $slice = $sa->fetch_by_region("chromosome", "19", "9533632", "9534232");

print $slice->seq;
