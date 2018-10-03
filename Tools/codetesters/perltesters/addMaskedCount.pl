#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Text::CSV;

my $spp = "mouse";
my $filename = "/home/pzs/histone/compositeprofile/trunk/datafiles/processed-ES-eistructures.csv";

Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ensembldb.ensembl.org', 
					      -user => 'anonymous');

my $sa = Bio::EnsEMBL::Registry->get_adaptor($spp, "core", "Slice");

my $delimiter = "\t";
# csv processing tool
my $csv = Text::CSV->new({sep_char => $delimiter, eol => "\n"});
my $coord_index = 4;


open (CSV, "<", $filename) or die $!;
# for each gene in file
while (<CSV>) {
	my @columns;
	my $coord;
	if ($csv->parse($_)) {
		@columns = $csv->fields();
		$coord = $columns[$coord_index];
	} else {
		my $err = $csv->error_input;
		print "Failed to parse line: $err";
	}
	$coord =~ /chr(\w+):(\d+)-(\d+)/;
	my $chr = $1;
	my $start = $2;
	my $end = $3;
	
	my $slice = $sa->fetch_by_region( 'chromosome', $chr, $start, $end );
	my $repeat_masked_slice = $slice->get_repeatmasked_seq();
	my $dna = $repeat_masked_slice->seq();
	my $basecount = length($dna);
	my $ncount = ($dna =~ s/N/$1/g);
	if(!$ncount) 
	{
		$ncount = 0;
	}
	if($basecount != ($end - $start) + 1)
	{
		die("base count mismatch! " . $chr . " " . $start . " " . $end . "\n");
	}
	my $codingcount = $basecount - $ncount;
# 	print "gene: " . $columns[0] . "\n";
# 	print "ei id: " . $columns[3] . "\n";
# 	print "ei: " . $columns[5] . "\n";
# 	print "coord: " . $coord . "\n";
# 	print "ncount: " . $ncount . "\n";
# 	print "regular count: " . $basecount . "\n";
# 	print "(end - start) + 1: ";
# 	print (($end - $start) + 1);
# 	print "\n";
# 	print "coding count: " . $codingcount . "\n";
# 
# 	print "\n";	
	push @columns, $codingcount;
	print join("\t", @columns) . "\n";
}
