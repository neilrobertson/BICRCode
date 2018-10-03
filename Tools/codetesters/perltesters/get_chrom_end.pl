#!/usr/bin/perl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";

$reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');

my $sa = $reg->get_adaptor('human', 'core', 'Slice');

my $basedir = "/home/pzs/genebuilds/human";
my $versionfile = $basedir . "/version.txt";
open VOUT, "> $versionfile" || die;

print VOUT $sa->dbc()->dbname() , "\n";
close VOUT;




my @chromosomes = (	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
					"12", "13", "14", "15", "16", "17", "18", "19", "20", "21", 
					"22", "X", "Y" );

print join("\t", ("chrom", "lastcoord")) . "\n";
foreach my $chr (@chromosomes) {           # make an array of chromosome names
#	print "working at chromosome: ". $chr . "\n";
	my $slice = $sa->fetch_by_region('chromosome',$chr) or die "Failed to slice chromosome $chr";
	
	print join("\t", ($chr, $slice->end)) . "\n";
	
}



