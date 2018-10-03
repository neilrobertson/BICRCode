#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;

my $spp = "human";
my $chr = 19;
my $start = 4000000;
my $end =   4005000;

Bio::EnsEMBL::Registry->load_registry_from_db(-host =>'ensembldb.ensembl.org', 
					      -user => 'anonymous');

my $sa = Bio::EnsEMBL::Registry->get_adaptor($spp, "core", "Slice");

my $slice = $sa->fetch_by_region( 'chromosome', $chr, $start, $end);

# method 1
my $repeat_masked_slice = $slice->get_repeatmasked_seq();
my $dna = $repeat_masked_slice->seq();
print "size of masked region according to get_repeatmasked_seq(): ";
print ($dna =~ s/N/$1/g);
print "\n";

# method 2
my @repeat_elements = @{$slice->get_all_RepeatFeatures()};
my @repeat_regions = ();
foreach my $re (@repeat_elements)
{
	my $feature = $re->transform("chromosome");
	my $restart = $feature->start;
	my $reend = $feature->end;
	my $name = $feature->display_id;
	my @thisregion = ($restart, $reend);
	push @repeat_regions, [ @thisregion ];
}

# my $first = @repeat_regions->[0];	
# $first->[0] = 1;
# 
print "before:\n";
foreach my $region (@repeat_regions)
{
	print $region->[0] . "\t" . $region->[1] . "\n";
}
print "\n";

my @merged_regions = ();
my $lastregion = 0;
foreach my $regionref (@repeat_regions)
{
	print "working at region (" . $regionref->[0] . ", " . $regionref->[1] . ")\n";
	# either we're at the start, or this region is entirely outside the last one
	if($lastregion == 0 or $regionref->[0] >= $lastregion->[1])
	{
		print "adding to list\n";
		push @merged_regions, $regionref ;
		$lastregion = $regionref;
	}
	else
	{
		# new region entirely inside the last one - discard
		if($regionref->[1] <= $lastregion->[1])
		{
			print "subsumed: skip\n\n";
			next;
		}
		else
		{
			$lastregion->[1] = $regionref->[1];
			print "updating last region - now (" . $lastregion->[0] . ", " . $lastregion->[1] . ")\n";
		}
	}
	print "\n";
}

if($merged_regions[0]->[0] < $start)
{
	$merged_regions[0]->[0] = $start;
}
	
my $regionsize = 0;
print "after\n";
foreach my $region (@merged_regions)
{
	print $region->[0] . "\t" . $region->[1] . "\n";
	$regionsize += ($region->[1] - $region->[0]) + 1;
}
print "total region size: " . $regionsize . "\n";
print "\n";

my $i;
my $repeatcount = 0;
for($i=$start; $i<=$end+1; $i++)
{
	my $found = 0;
	foreach my $re (@repeat_regions)
	{
		my $restart = $re->[0];
		my $reend = $re->[1];
		if($i<=$reend and $i >= $restart)
		{
			$found = 1;
			last;
		}
	}
	if($found == 1)
	{
		$repeatcount++;
	}
}

print "size of masked region with get_all_RepeatFeatures(): ";
print $repeatcount;
print "\n";

# my $i;
# my $repeatcount = 0;
# for($i=$start; $i<=$end; $i++)
# {
# 	my $found = 0;
# 	foreach my $re (@repeat_elements)
# 	{
# 		my $feature = $re->transform("chromosome");
# 		my $restart = $feature->start;
# 		my $reend = $feature->end;
# 		if($i<=$reend and $i >= $restart)
# 		{
# 			$found = 1;
# 			last;
# 		}
# 	}
# 	if($found == 1)
# 	{
# 		$repeatcount++;
# 	}
# }
# 
# print "size of masked region with get_all_RepeatFeatures(): ";
# print $repeatcount;
# print "\n";
