#!/usr/bin/perl

use strict;
use warnings;
use lib ("/home/pzs/src/ensembl-36/ensembl/modules/");
use Bio::SeqIO;
use Bio::EnsEMBL::Registry;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
											-user => 'anonymous',
											-dbname => 'homo_sapiens_core_36_35i');

my $ga = $db->get_GeneAdaptor();

my $gene = $ga->fetch_by_stable_id("ENSG00000176907");

print $gene->stable_id() . " " . $gene->start() . " " . $gene->end() . "\n";
