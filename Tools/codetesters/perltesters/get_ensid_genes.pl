#!/usr/bin/env perl

use Text::CSV;


my $filename = "/home/pzs/workspace/scriptsandbox/output/ensid-activations-rmapartitions2.csv";

# csv processing tool
my $csv = Text::CSV->new();
my $gene_index = 0;

my $geneout = "ensidgenes-rmapartitions2.csv";
my $exonout = "ensidexons-rmapartitions2.csv";

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;




## connect to database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');

my $gene_adaptor = $db->get_GeneAdaptor();
my %seen = ();

open(GENEOUT, ">$geneout") or die "could not open file ", $geneout;
open(EXONOUT, ">$exonout") or die "could not open file ", $exonout;

my $notfound = 0;

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
	my $gene = $gene_adaptor->fetch_by_stable_id($genename);
#	print "found " . scalar(@genes) . "\n";
	if($gene)
	{
		print "found!\n";

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
	else
	{
		print "not found!\n";
		$notfound = 1;
	}
}

if($notfound)
{
	print "some not found!";
}
