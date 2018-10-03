#!/usr/bin/perl -w

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);


my @chromosomes = (	"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
					"12", "13", "14", "15", "16", "17", "18", "19", "20", "21", 
					"22", "X", "Y" );


foreach(@chromosomes) { $sitecounts{$_} = 0; }
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

my $coord_system_name = "chromosome";
my $chr = "X";

print "working at chromosome: ". $chr . "\n";
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
	my $upstream = substr($chr_seq, $site+4, 100);          # $site starts at the 'G' in 'GTAC', so begin capture after 'C'
	my $downstream = substr($chr_seq, $site-100, 100);
	push @output, "$chr\t$truesite\t$upstream\t$downstream\n";
}

open FPOUT, "> tester.txt" || die;

foreach(@output) { print FPOUT $_; }

while( my($key, $value) = each(%sitecounts) ) {
	print "$key => $value\n";
}
