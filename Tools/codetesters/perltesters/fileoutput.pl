#!/usr/bin/perl

@output = ("first\n", "second\n", "third\n");

open FP, "> tester.txt" || die;

foreach(@output) { print FP $_; } 

close FP;
