#!/usr/bin/env perl

use Text::CSV;

use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");

#my $filename = "/home/pzs/publicdata/MAS5_RMA/output/RMA_output.csv";
#my $filename = "/home/pzs/publicdata/mikkelsen07mouse/celfiles/ES/output/RMA_output.csv";
#my $filename = "tester.csv";
#my $filename = "/home/pzs/histone/HISTONE_DATA/data/microarray/K562_mean_expression_rma2.txt";
my $filename = "/home/pzs/publicdata/CD4-expression-schones/resting/output/RMA_output.csv";
my $delimiter = "\t";
# csv processing tool
my $csv = Text::CSV->new({sep_char => $delimiter});
my $gene_index = 0;
my $species = "human";
my $celltype = "CD4";
my $ncbi;
if($species eq "mouse")
{	
	$ncbi = "NCBIm";
}
elsif($species eq "human")
{
	$ncbi = "NCBI";
}


#use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");
#use lib ("/home/pzs/src/ensembl/ensembl/modules/");

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;




## connect to database
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => "ensembldb.ensembl.org",
	-user => "anonymous"
);


my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );
my %seen = ();

my $fullversion = $gene_adaptor->dbc()->dbname;
my @versionlist = split(/_/, $fullversion);
$versioninfo = $ncbi . $versionlist[-1];
my $geneout = "genes-" . $species . "-" . $celltype . "-" . $versioninfo . ".csv";

print "writing into file: " . $geneout . "\n";

open(GENEOUT, ">$geneout") or die "could not open file ", $geneout;
#open(EXONOUT, ">$exonout") or die "could not open file ", $exonout;


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

	if(!$genename)
	{
		next;
	}
	print "looking for " . $genename . " ... ";
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
				my @line = ($genename, $gid, $gcoord, $strand, $transcriptname, $tcoord, $exonname, $ecoord);
				print GENEOUT join("\t", @line, "\n");
			}
		}

		$seen{$gid} = 1;
			
	}
	else
	{
		print "not found\n";
	}
}

print "size of seen hash ", scalar(keys(%seen)), "\n";
