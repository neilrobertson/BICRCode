#!/usr/bin/env perl

use strict;
use warnings;

use lib ("/home/pzs/src/ensembl-38/ensembl/modules/");
use Bio::SeqIO;

my $spp = "Mus Musculus";
use Bio::EnsEMBL::Registry;

my $reg = "Bio::EnsEMBL::Registry";

$reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');

my $sa = $reg->get_adaptor($spp, 'core', 'Slice');

print $sa->dbc()->dbname;
