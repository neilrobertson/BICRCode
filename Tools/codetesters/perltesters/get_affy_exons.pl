#!/usr/bin/env perl

use Text::CSV;

use lib ("/home/pzs/src/ensembl-55/ensembl/modules/");

my $filename = "/home/pzs/expressionarrays/kelly/output/RMA_output.csv";

my $delimiter = "\t";
my $species = "mouse";

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
my $exon_adaptor = $registry->get_adaptor( $species, 'Core', 'Exon' );
my $db = $gene_adaptor->db;
my $afa = $db->get_AffyFeatureAdaptor();


# open specified file
open (CSV, "<", $filename) or die $!;
# for each gene in file
while (<CSV>) {
	my $probename = "";
	if ($csv->parse($_)) {
		my @columns = $csv->fields();
		$probename = $columns[$gene_index];
	} else {
		my $err = $csv->error_input;
		print "Failed to parse line: $err";
	}

	if(!$probename)
	{
		next;
	}
	print "looking for " . $probename . " ... ";
 	my @hits = @{$afa->fetch_all_by_probeset($probename)};
 	print "found " . scalar(@hits) . "\n";

	foreach $hit (@hits)
	{
		@affyfeatures = @{$afa->fetch_all_by_AffyProbe($hit)};
		foreach $afeature (@affyfeatures) 
		{
			$feature = $afeature->transform('chromosome');
			print join("\t", ($feature->seq_region_name, $feature->start, $feature->end)) . "\n"
		}
	}

}
