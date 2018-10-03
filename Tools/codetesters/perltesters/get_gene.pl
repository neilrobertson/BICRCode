#!/usr/bin/perl

use Bio::EnsEMBL::Registry;

sub fetch_by_gene_id {
	my ($stable_id, $parameters) = @_;

	my ($host, $user) = ($parameters->{host}, $parameters->{user});

	my $registry = 'Bio::EnsEMBL::Registry';

	$registry->load_registry_from_db(
		-host => $host,
		-user => $user
	);
	
	my $gene_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );
	my $gene = $gene_adaptor->fetch_by_display_label($stable_id);
	
	return $gene;
	
}

my $parameters;

$parameters->{"host"} = "ensembldb.ensembl.org";
$parameters->{"user"} = "anonymous";

# scalar, because it's a reference
my $attributes;

$attributes = fetch_by_gene_id("IL13", $parameters)->get_all_Attributes();

foreach $attr (@$attributes)
{
	print $attr->name();
}
