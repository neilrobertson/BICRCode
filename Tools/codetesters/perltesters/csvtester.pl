#!/usr/bin/perl

use Text::CSV;

my $filename = "/home/pzs/histone/compositeprofile/trunk/datafiles/processed-eistructures.csv";

my $delimiter = "\t";
# csv processing tool
my $csv = Text::CSV->new({sep_char => $delimiter, eol => "\n"});
my $gene_index = 0;


open (CSV, "<", $filename) or die $!;
# for each gene in file
while (<CSV>) {
	my $genename = "";
	if ($csv->parse($_)) {
		my @columns = $csv->fields();
		$genename = $columns[$gene_index];
	} else {
		my $err = $csv->error_input;
		print "Failed to parse line: $err";
	}
	print "found genename: " . $genename . "\n";
}
