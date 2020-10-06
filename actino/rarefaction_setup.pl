#!/bin/env perl

use strict;
use warnings;

open my $ffh, '<', 'fullnetcc.tsv' or die $!;
while(<$ffh>){
    chomp;
    my ($cnum, $n) = split(/\t/, $_);
    foreach my $i (1..$n){
	print 'comp'.$cnum."\n"
    }
}
close $ffh;
my $singstart = 1000;
foreach (1..273){
    my $n = $singstart + $_;
    print 'comp'.$n."\n";
}
