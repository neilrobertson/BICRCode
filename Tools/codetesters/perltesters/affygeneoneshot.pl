#!/usr/bin/env perl

use Text::CSV;


use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;




## connect to database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');

my $gene_adaptor = $db->get_GeneAdaptor();

$genename = "117_at";

print "looking for " . $genename . "...";
my @genes = @{$gene_adaptor->fetch_all_by_external_name($genename)};
#	print "found " . scalar(@genes) . "\n";
if(scalar(@genes) == 1)
{
	print "found!\n";
	my $gene = $genes[0];
	my $gid = $gene->display_id();
	my $chrom = $gene->seq_region_name();
	my $gcoord = $chrom . ":" . $gene->start . "-" . $gene->end;
	my $strand = $gene->strand;
	my @geneline = ($genename, $gcoord, $strand);
	print join("\t", @geneline, "\n");

# 	my @transcripts = @{$gene->get_all_Transcripts()};
# 	# for each transcript
# 	foreach my $transcript (@transcripts) {
# 		my $transcriptname = $transcript->stable_id();
# 		# get exons
# 		my @exons = @{$transcript->get_all_Exons()};
# 		# dump lines to output file
# 		foreach my $exon (@exons) {
# 			my $exonname = $exon->stable_id();
# 			my $tcoord = $chrom . ":" . $transcript->start . "-" . $transcript->end;
# 			my $ecoord = $chrom . ":" . $exon->start . "-" . $exon->end;
# 			my @line = ($genename, $geneid, $gcoord, $strand, $transcriptname, $tcoord, $exonname, $ecoord);
# 			my @exonline = ($genename, $transcriptname, $tcoord, $exonname, $ecoord);
# 			print EXONOUT join("\t", @exonline, "\n");
# 		}
# 	}

}
else
{
	print "not found\n";
}
