#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);

my $spp = "human";
my $chr = $ARGV[0];
my $start = $ARGV[1];
my $end =   $ARGV[2];

Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ensembldb.ensembl.org', 
#Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ens-staging',
					      -user => 'anonymous');

my $sa = Bio::EnsEMBL::Registry->get_adaptor($spp, "core", "Slice");

my $slice = $sa->fetch_by_region( 'chromosome', $chr, $start, $end);

my $repeat_masked_slice = $slice->get_repeatmasked_seq();

# get repeat masked sequence:
my $dna = $repeat_masked_slice->seq();
#$dna = $repeat_masked_slice->subseq( 1, 1000 );

print $dna . "\n";
