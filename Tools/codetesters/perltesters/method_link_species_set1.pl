#
# This script counts the number of MethodLinkSpeciesSets stored in
# the compara database
#
use strict;

use Bio::EnsEMBL::Registry;

# Auto-configure the registry
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host=>"ensembldb.ensembl.org", -user=>"anonymous");

# Get the Compara Adaptor for MethodLinkSpeciesSets
my $mlssa = Bio::EnsEMBL::Registry->get_adaptor(
    "Multi", "compara", "MethodLinkSpeciesSet");

# fetch_all() method returns a ref. to an array
my $all_mlss = $mlssa->fetch_all();

# @$all_mlss returns the array referenced by $all_mlss. The
# scalar value of the array is the number of elements in the
# array, i.e. the number of MethodLinkSpeciesSets in the
# database
print scalar(@$all_mlss),"\n";
