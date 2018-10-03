#!/usr/bin/perl

use strict;

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);

my $dba = $registry->get_DBAdaptor("human", "core");

print $dba->dbc()->dbname() . "\n";

