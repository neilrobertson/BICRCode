#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";

$reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');

my $sa = $reg->get_adaptor('human', 'core', 'Slice');
print $sa->dbc()->dbname() , "\n";

#my $slice = $sa->fetch_by_region('chromosome', '20', 1, 10000000);

#print $slice->get_seq_region_id(), "\n";

# my $output = Bio::SeqIO->new( -file=>'>2a.txt', -format=>'FASTA');
# $output->write_seq($slice);


