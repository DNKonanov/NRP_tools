#!/bin/env perl

use strict;
use warnings;

my $tot = 0;
while(<>){
    chomp;
    my @a = split(/\t/, $_);
    $tot += $a[2];
}
print "$tot\n";
