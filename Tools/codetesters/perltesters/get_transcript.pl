#!/usr/bin/perl

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
	
	return $transcript->seq()->seq();
	
}

my $parameters;

$parameters->{"host"} = "ensembldb.ensembl.org";
$parameters->{"user"} = "anonymous";


print fetch_by_transcript_id("ENST00000368346", $parameters);
