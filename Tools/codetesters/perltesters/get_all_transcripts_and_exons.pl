#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";

$reg->load_registry_from_db( -host => 'ensembldb.ensembl.org', -user => 'anonymous');

my $species = "human";

my $sa = $reg->get_adaptor($species, 'core', 'Slice');



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
		my $strand = $gene->strand();
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
				my @line = ($genename, $gcoord, $strand, $transcriptname, $tcoord, $exonname, $ecoord);
				print  join("\t", @line, "\n");
			}
		}
	}	
}



