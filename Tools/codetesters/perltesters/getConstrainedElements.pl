#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::AlignIO;

#
# Simple example to retrieve the constrained elements from the PECAN 10 way or
# ORTHEUS 23 way alignment
#

my $spp = "Mus musculus";
my $chr = $ARGV[0];
my $start = $ARGV[1];
my $end =   $ARGV[2];


#EXAMPLE FOR ENSEMBL 53
#****CONSTRAINED ELEMENTS HAVE BEEN MOVED FROM THE GENOMIC_ALIGN AND GENOMIC_ALIGN_BLOCK TABLES IN A CONSTRAINED_ELEMENT TABLE****

#Species in EPO 31 way alignment
my $spp_list = ["Rattus norvegicus","Homo sapiens","Macaca mulatta","Echinops telfairi","Oryctolagus cuniculus","Pan troglodytes","Canis familiaris","Tupaia belangeri","Erinaceus europaeus","Otolemur garnettii","Spermophilus tridecemlineatus","Myotis lucifugus","Sorex araneus","Mus musculus","Microcebus murinus","Pongo pygmaeus","Equus caballus","Bos taurus","Felis catus","Ochotona princeps","Cavia porcellus","Gorilla gorilla","Choloepus hoffmanni","Procavia capensis","Tursiops truncatus","Loxodonta africana","Tarsius syrichta","Dipodomys ordii","Vicugna pacos","Pteropus vampyrus","Dasypus novemcinctus"];

my $version = 54;
my $original_alignment_mlss_id = 364; #mlssid from alignments from which the constrained elements were derived.

Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ensembldb.ensembl.org', 
#Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ens-staging',
					      -user => 'anonymous');


#Create slice from $spp, $chr, $start and $end
my $query_slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($spp, "core", "Slice");
my $query_slice = $query_slice_adaptor->fetch_by_region("chromosome",$chr, $start, $end);

# Getting the MethodLinkSpeciesSet adaptor: 
my $method_link_species_set_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet');

#Get constrained element method_list_species_set
my $ce_mlss = $method_link_species_set_adaptor->fetch_by_method_link_type_registry_aliases("GERP_CONSTRAINED_ELEMENT", $spp_list);

#Get the mlss object from the alignement from which the constrained elements were generated (31-way ORTHEUS)
my $orig_mlss = $method_link_species_set_adaptor->fetch_by_dbID( $original_alignment_mlss_id );


#Get constrained_element adaptor
my $ce_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'ConstrainedElement');

#Fetch all constrained elements
my $cons = $ce_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($ce_mlss,$query_slice);

print "Number of constrained elements: " . @$cons . "\n";

#Print out information
#Note: where constrained elements occur in overlapping genomic_align_blocks there will be ambiguities
#in aassociating an alignment with the correct constrined_element_id. 
foreach my $ce (@$cons) {
    print "dbID:" . $ce->dbID . " from:" . ($ce->slice->start + $ce->start - 1 ) . " to:" . 
	($ce->slice->start + $ce->end - 1) . "\n";
}

