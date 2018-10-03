#!/usr/bin/perl

use strict;
use warnings;
use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");
use Bio::EnsEMBL::Registry;

use Bio::SeqIO;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');

my $sa = $db->get_SliceAdaptor();



my @chromosomes = (	
					"1", 
 					"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
 					"12", "13", "14", "15", "16", "17", "18", "19",  
 					"20", "21", "22",
					"X", "Y" 
					);

foreach my $chr (@chromosomes) {           # make an array of chromosome names
	my $slice = $sa->fetch_by_region('chromosome',$chr) or die "Failed to slice chromosome $chr";
	
	my @genes = @{ $slice->get_all_Genes() };
	foreach my $gene (@genes) {
		my $chrom = $gene->seq_region_name();
		my $gcoord = $chrom . ":" . $gene->start . "-" . $gene->end;
		my $genename = $gene->stable_id();
		# get transcripts
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
				my @line = ($genename, $gcoord, $transcriptname, $tcoord, $exonname, $ecoord);
				print  join("\t", @line, "\n");
			}
		}
	}	
}



