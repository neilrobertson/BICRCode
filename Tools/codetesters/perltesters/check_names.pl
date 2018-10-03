#!/usr/bin/env perl

use Text::CSV;


my $filename = "/home/pzs/workspace/scriptsandbox/output/affyactivations.csv";
# csv processing tool
my $csv = Text::CSV->new();
my $gene_index = 0;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;




## connect to database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');

my $gene_adaptor = $db->get_GeneAdaptor();
my %seen = ();

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

	print "looking for " . $genename . "...";
	my @genes = @{$gene_adaptor->fetch_all_by_external_name($genename)};
#	print "found " . scalar(@genes) . "\n";
	if(scalar(@genes) == 1)
	{
		print "found!\n";
		my $gene = $genes[0];
		my $gid = $gene->display_id();
		if(!$seen{$gid})
		{
			$seen{$gid} = 1;
		}
	}
	else
	{
		print "not found\n";
	}
}

print "size of seen hash ", scalar(keys(%seen)), "\n";
