#!/bin/env perl

use strict;
use warnings;

my %tograb = ();
open my $pfh, '<', 'putclust.tsv' or die $!;
while(<$pfh>){
    chomp;
    my($contig, $k, $doms) = split(/\t/, $_);
    $tograb{$k} = 1;
}
close $pfh;
open my $ofh, '<', 'all.orf.tsv' or die $!;
while(<$ofh>){
    chomp;
    my ($contig, $start, $end, $spec, $k, $orf) = split(/\t/, $_);
    if(exists $tograb{$k}){
	print join("\t", $contig, $start, $end, $spec, $k, $orf)."\n";
    }
}
close $ofh;
