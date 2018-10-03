#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ensembldb.ensembl.org', 
#Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ens-staging',
					      -user => 'anonymous');

my $spp = "mouse";

my $ga = Bio::EnsEMBL::Registry->get_adaptor($spp, "core", "Gene");
my $sa = Bio::EnsEMBL::Registry->get_adaptor($spp, "core", "Slice");

sub dumpStopCodonsOneGeneID
{
	my $geneid = shift;
	my $gene = $ga->fetch_by_stable_id( $geneid );
	my $chrom = $gene->seq_region_name();
	my $strand = $gene->strand;
	
	my $extname = $gene->external_name;

	my @transcripts = @{$gene->get_all_Transcripts()};

	foreach my $transcript (@transcripts) {
		my $transcriptname = $transcript->stable_id();
		# length of all exons - different to genome coverage...
		my $tlength = $transcript->length();
		# get exons
		my @exons = @{$transcript->get_all_Exons()};
		my $lastexon = $exons[-1];
		my $lestart = $lastexon->start;
		my $leseq = $lastexon->seq();

		my $pos;
	
		my $scstart;
		my $scend;
	
		if($strand == -1)
		{
			$scstart = $transcript->coding_region_start;
			$scend = $scstart + 2;
		}
		else
		{
			$scend = $transcript->coding_region_end;
			$scstart = $scend - 2;
		}

		my $slice = $sa->fetch_by_region( "chromosome", $chrom, $scstart, $scend );
		my $stopcodon;
		if($strand == -1)
		{
			my $scsec = Bio::Seq->new(-seq => $stopcodon = $slice->seq, alphabet => 'dna');
			$stopcodon = $scsec->revcom->seq;
		}
		else
		{
			$stopcodon = $slice->seq;
		}
		
 		my $scrange = "chr" . $chrom . ":" . $scstart . "-" . $scend;
 		my $insertion;
		if($strand == -1)
		{
			$insertion = "chr" . $chrom . ":" . join("-", $scend, $scend+1);
		}
		else
		{
			$insertion = "chr" . $chrom . ":" . join("-", $scstart-1, $scstart);
		}
		# dump lines to output file
		my @line = ($extname, $geneid, $strand, $transcriptname, $tlength, $stopcodon, $scrange, $insertion);
		print join("\t", @line, "\n");
	}
}

sub dumpStopCodonsOneGeneName
{
	my $genename = shift;
	
	my @genes = $ga->fetch_all_by_external_name($genename);
	
	if(scalar(@genes) == 1)
	{
		my $geneid = $genes[0];
		dumpStopCodonsOneGeneID($geneid);
	}
	else
	{
		print "more than one id for genename: " . $genename . "\n";
	}
}

my @line = ("gene name", "gene id", "strand", "transcript id", "transcript length", "stop codon", "stop codon range", "insertion position");
print join("\t", @line, "\n");

#my $geneid = $ARGV[0];
my @geneids = (	"ENSMUSG00000044791", "ENSMUSG00000021488", "ENSMUSG00000033326",
				"ENSMUSG00000054611", "ENSMUSG00000029475", "ENSMUSG00000061589");
#my @geneids = (	"ENSMUSG00000061589" );

if(not @ARGV)
{
	@ARGV = @geneids;
}

for (@ARGV)
{
	dumpStopCodonsOneGeneID($_);
}
