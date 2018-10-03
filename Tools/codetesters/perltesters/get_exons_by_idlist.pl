#!/usr/bin/env perl


use Text::CSV;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");

# csv processing tool
my $csv = Text::CSV->new();
my $gene_index = 0;

my $filename = "/home/pzs/histone/IAN_newbuild.csv";
#my $filename = "/home/pzs/histone/HISTONE_DATA/ENCODE_IDSu.txt";
#my $filename = "/home/pzs/codetesters/perltesters/missing.txt";

# use Bio::EnsEMBL::Registry;
# my $reg = "Bio::EnsEMBL::Registry";
# $reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');


## connect to database
use Bio::EnsEMBL::DBSQL::DBAdaptor;
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');
my $vdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(	-host => 'ensembldb.ensembl.org',
												-user => 'anonymous',
												-dbname => 'homo_sapiens_vega_36_35i');
my $ga = $db->get_GeneAdaptor();
my $va = $vdb->get_GeneAdaptor();

# my $ga = $reg->get_adaptor('human', 'core', 'Gene');
# my $va = $reg->get_adaptor('human', 'vega', 'Gene');

open(IN, $filename);

# for each gene in file
while (<IN>) {
 	my $genename;
	if ($csv->parse($_)) {
		my @columns = $csv->fields();
		$genename = $columns[$gene_index];
	} else {
		my $err = $csv->error_input;
		print "Failed to parse line: $err";
	}

# 	chomp;
# 	my $genename = $_;
# 	next if($genename eq "");
# 	my $gene="";
# 
#	print "looking for ", $genename, "...\n" ;

	

# 	if(substr($genename, 0, 3) eq "ENS")
# 	{
# 		$gene = $ga->fetch_by_stable_id($genename);
# 	}
# 	else
# 	{
# 		$gene = $va->fetch_by_stable_id($genename);
# 	}
# 
	my @coregenes = @{$ga->fetch_all_by_external_name($genename)};
	if(scalar(@coregenes > 0)) 
	{ 
#		print "found in core\n";
		$gene = $coregenes[0];
	}
	else 
	{
		my @vegagenes = @{$va->fetch_all_by_external_name($genename)};
		if(scalar(@vegagenes > 0)) 
		{ 
#			print "found in vega\n";
			$gene = $vegagenes[0];
		}
		else
		{
			print "not found: ", $genename, "\n";
			next;
		}
	}	
#	print "gene: ", $genename, " found (core, vega) (", scalar(@coregenes), ",", scalar(@vegagenes), ")", "\n";
	
	my $geneextname = $gene->external_name();
	# get transcripts
	my @transcripts = @{$gene->get_all_Transcripts()};
	# for each transcript
	foreach my $transcript (@transcripts) {
		my $transcriptname = $transcript->stable_id();
		my $chrom = $transcript->seq_region_name();
		# get exons
		my @exons = @{$transcript->get_all_Exons()};
		# dump lines to output file
		foreach my $exon (@exons) {
			my $exonname = $exon->stable_id();
			my $tcoord = $chrom . ":" . $transcript->start . "-" . $transcript->end;
			my $ecoord = $chrom . ":" . $exon->start . "-" . $exon->end;
			my @line = ($geneextname, $transcriptname, $tcoord, $exonname, $ecoord);
			print  join("\t", @line, "\n");
		}
	}

}

