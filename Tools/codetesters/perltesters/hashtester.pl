#!/usr/bin/perl

my %seen = ();

$seen{'key1'} = 1;
$seen{'key2'} = 1;

if(!$seen{'key3'})
{
	print "yes!\n";
}

print "size of hash ", scalar(keys(%seen)), "\n";
