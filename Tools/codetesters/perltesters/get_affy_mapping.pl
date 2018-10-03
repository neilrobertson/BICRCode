#!/usr/bin/env perl

use Text::CSV;

use lib ("/home/pzs/src/ensembl-54/ensembl/modules/");

my $geneout = "csvoutput/complete-affy-mapping-HG_U133_A-NCBI36.csv";

my $delimiter = "\t";
my $filename = "/home/pzs/expressionarrays/K562/output/RMA_output.csv";
my $csv = Text::CSV->new({sep_char => $delimiter});
my $species = "human";

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;




## connect to database
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);


my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );

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
#	print "found " . scalar(@genes) . "\n";
	if(scalar(@genes) == 1)
	{
		my $gene = $genes[0];
		my $gid = $gene->display_id();
		print "$gid\n";
		print GENEOUT join(",", ($genename, $gid)) . "\n";
	}
	else
	{
		print "not found\n";
	}
}
