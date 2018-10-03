#!/software/bin/perl -w


use strict;
use warnings;

use lib ("/software/microarray/lib/ensembl/bioperl-1.2.3/");
use lib ("/software/microarray/lib/ensembl/ensembl-core-36/modules/");


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;


## genes of interest are those with more than 3 exons and greater than 6000bp in length
my ($max_no_exons, $max_length) = (3, 6000);


## connect to database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'ensembldb.ensembl.org',
																						-user => 'anonymous',
																						-dbname => 'homo_sapiens_core_36_35i');


## get a slice adaptor
my $slice_adaptor = $db->get_SliceAdaptor();


## foreach chromosome in turn
foreach my $chr (1..22, "X", "Y", "M") {

	## get the chromosome object
	my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr);
	

	## for all genes on the chromosome
	foreach my $gene (@{$slice->get_all_Genes_by_type('protein_coding', 'ensembl')}) {


		## if gene is long enough
		if (($gene->end - $gene->start) >= $max_length) {

			my %occurrence_of_exon;

			## count the occurrence of each exon in all the transcripts
			## NB: exons occuring in all transcripts are canonical

			foreach my $transcript (@{$gene->get_all_Transcripts}) {

				foreach my $exon (@{$transcript->get_all_Exons}) {

					if (exists $occurrence_of_exon{$exon->stable_id}) {
						$occurrence_of_exon{$exon->stable_id} += 1;
					} else {
						$occurrence_of_exon{$exon->stable_id} = 1;
					}
				}
			}

			print $gene->stable_id."\n";


			foreach my $exon (@{$gene->get_all_Exons}) {
				print $exon->stable_id." ".$exon->start." ".$exon->end." ".$exon->strand." ".$occurrence_of_exon{$exon->stable_id}."\n";
			}

			print "\n";
		}
	}
}

close FPOUT;
