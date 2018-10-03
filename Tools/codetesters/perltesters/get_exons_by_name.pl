#!/usr/bin/perl

use strict;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");

use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');

my $ga = $db->get_GeneAdaptor();

my $filename = "/home/pzs/codetesters/pytesters/exongeneids.txt";

open(FH, $filename);
my @genenames = ();

while(<FH>)
{
	chomp($_);
	if($_){ push(@genenames, $_); }
}

close(FH);

foreach my $genename (@genenames)
{
#	print "working at gene: ", $genename, "\n";
	my $gene = $ga->fetch_by_stable_id($genename);
	if($gene)
	{
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
	else
	{
		print STDERR "could not find: ", $genename, "\n";
	}
}

