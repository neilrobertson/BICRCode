#!/usr/bin/perl

use Bio::EnsEMBL::Registry;

my $parameters;

$parameters->{"host"} = "ensembldb.ensembl.org";
$parameters->{"user"} = "anonymous";

my ($host, $user) = ($parameters->{host}, $parameters->{user});

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => $host,
	-user => $user
);

my $gene_adaptor = $registry->get_adaptor( 'Mouse', 'Core', 'Gene' );

my $file = "/home/pzs/histone/all_genes.csv";
my $outfile = "all_mouse_genes.csv";

open(NEWCSV, ">$outfile");
open(MYINPUTFILE, "<$file");
while(<MYINPUTFILE>)
{
	my($line) = $_;
	chomp($line);
	$line =~ s/"//g;
	($gene_id,$region,$chrom,$start,$end,$strand) = split(",", $line);
	
	my $gene = $gene_adaptor->fetch_by_display_label($gene_id);
	if(not $gene)
	{
		print "no match for id: $gene_id\n";
		next;
	}
	my $mousestart = $gene->start;
	my $mouseend = $gene->end;

	print NEWCSV "$gene_id,$region,$chrom,$mousestart,$mouseend,$strand\n";
}

close(MYINPUTFILE);
close(NEWCSV);
