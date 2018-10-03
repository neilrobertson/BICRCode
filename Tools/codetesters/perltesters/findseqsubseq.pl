#!/usr/bin/perl

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");

use Bio::EnsEMBL::DBSQL::DBAdaptor;
## connect to database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
																						-user => 'anonymous',
																						-dbname => 'homo_sapiens_core_36_35i');


## get a slice adaptor
my $slice_adaptor = $db->get_SliceAdaptor();

my $chrom = "1";
my $lower = 47262288;
my $upper = 47518923;
my $sitecount = 0;

my $chr_slice = $slice_adaptor->fetch_by_region("chromosome", $chrom, $lower, $upper) or die "Failed to get slice";
my $chr_seq = $chr_slice->seq();
print length($chr_seq)." bases read.\n";        # just to check
my ($site, $prev) = (0, 0, -1);
until ($site == -1) {
 	$sitecount++;            # record number of sites per chromosome
 	$site = index($chr_seq, 'GTAC', $prev + 1);
 	last if ($site == -1);
 	$prev = $site;
 	my $truesite = $site + $lower;
	my $trueupper = $truesite + 4;
 	push @output, "chr$chrom",":$truesite","-","$trueupper\n";
}
print "site count: " . $sitecount . "\n";
 
open FPOUT, "> tester.txt" || die;

foreach(@output) { print FPOUT $_; }
close FPOUT;
