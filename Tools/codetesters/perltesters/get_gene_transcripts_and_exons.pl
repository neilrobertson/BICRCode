#!/usr/bin/perl -w

use Bio::EnsEMBL::Registry;
use Text::CSV;

my $filename = "/home/pzs/histone/newourmousegene.csv";

# csv processing tool
my $csv = Text::CSV->new();
my $gene_index = 4;

my $registry = 'Bio::EnsEMBL::Registry';

my $parameters;

$parameters->{"host"} = "ensembldb.ensembl.org";
$parameters->{"user"} = "anonymous";

my ($host, $user) = ($parameters->{host}, $parameters->{user});

$registry->load_registry_from_db(
	-host => $host,
	-user => $user
);

sub fetch_by_transcript_id {
	my ($stable_id) = @_;

	my $transcript_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Transcript' );
	my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);
	
	return $transcript;
}

sub fetch_by_gene_id {
	my ($stable_id) = @_;

	my $gene_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
	my $gene = $gene_adaptor->fetch_by_stable_id($stable_id);
	
	return $gene;
}


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
	my $gene = fetch_by_gene_id($genename);
	if (not($gene)) {
		print STDERR "could not find gene " . $genename . "\n";
		#die("could not find gene " . $genename . "\n");
		next;
	}
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
