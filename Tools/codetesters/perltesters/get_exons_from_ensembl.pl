#!/usr/bin/perl -w


use strict;
use warnings;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;




## connect to database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');


## get a slice adaptor
my $slice_adaptor = $db->get_SliceAdaptor();

# csv processing tool
my $csv = Text::CSV->new();
my $gene_index = 0;

## foreach chromosome in turn
foreach my $chr (1..22, "X", "Y", "M") {

	## get the chromosome object
	my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr);
	

	## for all genes on the chromosome
	foreach my $gene (@{$slice->get_all_Genes_by_type('protein_coding', 'ensembl')}) {

		my $genename = $gene->display_id();
		
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
				my @line = ($genename, $transcriptname, $tcoord, $exonname, $ecoord);
				print  join("\t", @line, "\n");
			}
		}
		


	}
}

