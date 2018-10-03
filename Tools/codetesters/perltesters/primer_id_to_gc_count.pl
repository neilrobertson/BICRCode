#!/usr/bin/env perl


use Text::CSV;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");

# csv processing tool
my $csv = Text::CSV->new();
my $gene_index = 0;

my $filename = "human_ids.txt";
#my $filename = "/home/pzs/histone/HISTONE_DATA/ENCODE_IDSu.txt";
#my $filename = "/home/pzs/codetesters/perltesters/missing.txt";

# use Bio::EnsEMBL::Registry;
# my $reg = "Bio::EnsEMBL::Registry";
# $reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');

$p_index = 0;
$coord_index = 1;

## connect to database
use Bio::EnsEMBL::DBSQL::DBAdaptor;
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');
my $sa = $db->get_SliceAdaptor();


open(IN, $filename);

# for each gene in file
while (<IN>) {
 	my $genename;
	if ($csv->parse($_)) {
		my @columns = $csv->fields();
		$pname = $columns[$p_index];
		$coord = $columns[$coord_index];
	} else {
		my $err = $csv->error_input;
		print "Failed to parse line: $err";
		next;
	}

	my @coordsplit1 = split(":", $coord);

	my $chrom = substr($coordsplit1[0], 3);

	my @coordsplit2 = split("-", $coordsplit1[1]);
	my $start = @coordsplit2[0];
	my $end = @coordsplit2[1];


	my $slice = $sa->fetch_by_region("chromosome", $chrom, $start, $end);
	my $seq = $slice->seq();
#	print $seq, "\n";
	
	# impossibly obscure idiom. God I hate Perl
	$gccount = ($seq =~ tr/GC//);
	$tacount = ($seq =~ tr/TA//);

	$prop = $gccount / $tacount;

	print $pname, "\t", $prop, "\n";
# 	my @line = ($pname, $prop);
# 	print join("\t", @line, "\n")
	
}

