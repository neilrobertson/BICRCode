#!/usr/bin/env perl

use strict;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::EnsEMBL::Registry;


my $registry = 'Bio::EnsEMBL::Registry';

my $dba = $registry->load_all();

# 
# $registry->load_registry_from_db(
# 	-host => "ensembldb.ensembl.org",
# 	-user => "anonymous"
# );
# my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
# 						    -user   => 'anonymous',
# 						    -dbname => 'homo_sapiens',
# 						    -host   => 'ensembldb.ensembl.org',
# 						    -driver => 'mysql' );

print $dba->dbc->dbname . "\n";

# my $ga = $db->get_adaptor("Homo Sapiens","core","Gene");
# 
# my $genename = $ARGV[0];
# 
# my @genes = @{$ga->fetch_all_by_external_name($genename)};
# if(scalar(@genes) == 1)
# {
# 	my $gene = $genes[0];
# 	print $genename, ",", $gene->strand, ",", "chr", $gene->seq_region_name, ":", $gene->start, "-", $gene->end, "\n"
# }
# 
