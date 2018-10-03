#!/usr/bin/env perl

use Text::CSV;

use lib ("/home/pzs/src/ensembl-54/ensembl/modules/");

my $filename = "siHRA_significant_genes.csv";
my $geneout = "siHRA_with_ensids.csv";
my $species = "human";

my $delimiter = ",";
# csv processing tool
my $csv = Text::CSV->new({sep_char => $delimiter});

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

## connect to database
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);


my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );
my %seen = ();

open(GENEOUT, ">$geneout") or die "could not open file ", $geneout;

# open specified file
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

	print "looking for " . $genename . " ... ";
	my @genes = @{$gene_adaptor->fetch_all_by_external_name($genename)};

	if(scalar(@genes) == 1)
	{
		my $gene = $genes[0];
		my $gid = $gene->display_id();
		print "$gid\n";
		my $newline = $gid . "," . $_;
		print GENEOUT $newline;
		$seen{$gid} = 1;
	}
	else
	{
		print "not found!\n";
	}
}
