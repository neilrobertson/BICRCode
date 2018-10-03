#!/usr/bin/perl

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);

$species = "mouse";

$slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );

my @chromosomes = (	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
					"12", "13", "14", "15", "16", "17", "18", "19", 
# 					"20", "21", 
# 					"22", 
					"X", "Y" );

foreach(@chromosomes) { $sitecounts{$_} = 0; }
my $coord_system_name = "chromosome";

my $markervalue = 100;

foreach $chr (@chromosomes) {           # make an array of chromosome names
	print "working at chromosome: ". $chr . "\n";
	my $fname = "gtacwiggles/" . $species . "/" . "chr" . $chr . ".txt";
	open FPOUT,  "> $fname" || die;
	print FPOUT "browser hide all\n";
	print FPOUT "browser full altGraph\n";
	print FPOUT 'track type=wiggle_0 name="GTAC: chr' . $chr .'" description="GTAC: chr' . $chr .  '"' . "visibility=full color=0,0,0 viewLimits=0:2 yLineMark=0 yLineOnOff=on useScore=0 priority=0 offset=0 autoscale=on gridDefault=off maxheightPixels=128:128:5 graphType=bar\n";

	print "dumping into file " . $fname . "\n";
	my $chr_slice = $slice_adaptor->fetch_by_region($coord_system_name,$chr) or die "Failed to slice chromosome $chr";
	my $chr_seq = $chr_slice->seq();
	print length($chr_seq)." bases read.\n";        # just to check
	my ($site, $prev) = (0, 0, -1);
	until ($site == -1) {
		$sitecounts{$chr}++;            # record number of sites per chromosome
		$site = index($chr_seq, 'GTAC', $prev + 1);
		last if ($site == -1);
		$prev = $site;
		my $truesite = $site + 1;       # because $chr_seq starts at 0, according to index()
		my $siteoffset = $truesite + 3;
		print FPOUT "chr$chr\t$truesite\t$siteoffset\t$markervalue\n";
	}
	print "site count: " . $sitecounts{$chr} . "\n";
}


foreach(@output) { print FPOUT $_; }
close FPOUT;

while( my($key, $value) = each(%sitecounts) ) {
	print "$key => $value\n";
}

