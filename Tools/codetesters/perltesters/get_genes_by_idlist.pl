#!/usr/bin/env perl

use Bio::EnsEMBL::Registry;
use Text::CSV;

#my $filename = "/home/pzs/histone/IAN_newbuild.csv";
my $filename = "/home/pzs/codetesters/pytesters/exongeneids.txt";

$parameters->{"host"} = "ensembldb.ensembl.org";
$parameters->{"user"} = "anonymous";

my ($host, $user) = ($parameters->{host}, $parameters->{user});

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => $host,
	-user => $user
);

my $gene_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );

open(FH, $filename);
my @genenames = ();

while(<FH>)
{
	chomp($_);
	push(@genenames, $_);
}

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

	print "looking for " . $genename . "...";
	my @genes = @{$gene_adaptor->fetch_by_stable_id($genename)};
	print "found " . scalar(@genes) . "\n";
}
