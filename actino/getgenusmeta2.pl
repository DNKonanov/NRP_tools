#!/bin/env perl

use strict;
use warnings;

my %gen = ();
open my $gfh, '<', 'genuscount.tsv' or die $!;
while(<$gfh>){
    chomp;
    my ($c, $g) = split(/\t/, $_);
    chomp(my $a = `ls ~/genomes/20160811-publicactinos/fna/*|grep $g |wc -l`);
    print join(',', $g, sprintf("%.3f", log10($a)))."\n";
}
close $gfh;

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}
