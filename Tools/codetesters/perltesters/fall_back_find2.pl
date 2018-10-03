#!/usr/bin/env perl


use Text::CSV;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");

# csv processing tool
my $csv = Text::CSV->new();
my $name_index = 0;
my $id_index = 1;

my $filename = "namesandids.csv";
#my $filename = "/home/pzs/histone/HISTONE_DATA/ENCODE_IDSu.txt";
#my $filename = "/home/pzs/codetesters/perltesters/missing.txt";
my $geneout = "biomartgenes.csv";
my $exonout = "biomartexons.csv";

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

open(IN, $filename) or die "could not open file ", $filename;
open(GENEOUT, ">$geneout") or die "could not open file ", $geneout;
open(EXONOUT, ">$exonout") or die "could not open file ", $exonout;

my $misscount = 0;

my @line = ("gene name", "geneid", "gene coord", "strand", "transcript name", "transcript coord", "exon name", "exon coord");
# for each gene in file
while (<IN>) {
 	my $genename;
	my $geneid;
	if ($csv->parse($_)) {
		my @columns = $csv->fields();
		$genename = $columns[$name_index];
		$geneid = $columns[$id_index];
	} else {
		my $err = $csv->error_input;
		print "Failed to parse line: $err";
	}

	chomp($geneid);
	chomp($genename);

 	if(substr($geneid, 0, 3) eq "ENS")
 	{
 		$gene = $ga->fetch_by_stable_id($geneid);
 	}
 	else
 	{
 		$gene = $va->fetch_by_stable_id($geneid);
 	}

	if($gene)
	{
#		print "found id ", $geneid, " in either vega or ensembl\n";
	}
	else
	{
		# look for name
		@genes = @{$ga->fetch_all_by_external_name($genename)};
		if(scalar(@genes > 0))
		{
			$gene = $genes[0];
#			print "found name ", $genename, " as external name\n";
		}
		else
		{
			my @line = ($genename, $geneid, "not found!");
			print join("\t", @line, "\n");
			$misscount++;
			next;
		}
	}

	my $chrom = $gene->seq_region_name();
	my $gcoord = $chrom . ":" . $gene->start . "-" . $gene->end;
	my $strand = $gene->strand;
	my @geneline = ($genename, $gcoord, $strand);
	print GENEOUT join("\t", @geneline, "\n");

	my @transcripts = @{$gene->get_all_Transcripts()};
	# for each transcript
	foreach my $transcript (@transcripts) {
		my $transcriptname = $transcript->stable_id();
		# get exons
		my @exons = @{$transcript->get_all_Exons()};
		# dump lines to output file
		foreach my $exon (@exons) {
			my $exonname = $exon->stable_id();
			my $tcoord = $chrom . ":" . $transcript->start . "-" . $transcript->end;
			my $ecoord = $chrom . ":" . $exon->start . "-" . $exon->end;
			my @line = ($genename, $geneid, $gcoord, $strand, $transcriptname, $tcoord, $exonname, $ecoord);
			my @exonline = ($genename, $transcriptname, $tcoord, $exonname, $ecoord);
			print EXONOUT join("\t", @exonline, "\n");
		}
	}

}

print "missed: ", $misscount, "\n";
