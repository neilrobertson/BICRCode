#!/usr/bin/perl -w

use strict;
use Bio::EnsEMBL::Registry;

sub fetch_by_transcript_id {
	my ($stable_id, $parameters) = @_;

	my ($host, $user) = ($parameters->{host}, $parameters->{user});

	my $registry = 'Bio::EnsEMBL::Registry';

	$registry->load_registry_from_db(
		-host => $host,
		-user => $user
	);
	
	my $transcript_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
	my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);
	
	return $transcript;
}

my $parameters;

$parameters->{"host"} = "ensembldb.ensembl.org";
$parameters->{"user"} = "anonymous";

my $transcript = fetch_by_transcript_id("ENST00000368346", $parameters);

my $chrom = $transcript->seq_region_name();

my @exons = @{$transcript->get_all_Exons()};
	

# print $#{$exons};
# print "\n";

foreach my $exon (@exons)
{
	print $chrom . ":" . $exon->start . "-" . $exon->end . "\n";
}
