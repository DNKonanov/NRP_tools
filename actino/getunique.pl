#!/bin/env perl

use strict;
use warnings;

my $n = 0;
while(<>){
    chomp;
    my ($count, $doms) = split(/\t/, $_);
    $doms =~ s/\|//g;
    if(length($doms) >= 3){
	$n += 1;
    }
}
print "$n\n";
