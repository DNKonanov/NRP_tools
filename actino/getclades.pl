#!/bin/env perl

use strict;
use warnings;

my %taxof = ();
chomp(my $q = `grep '>' ~/genomes/20160811-publicactinos/fna/*`);
foreach my $res (split(/\n/, $q)){
    my ($source, @rest) = split(/:/, $res);
    my @sd = split(/\//, $source);
    my @t = split(/_/, $sd[-1]);
    $rest[0] =~ s/^>(\S+).+$/$1/;
    $taxof{$rest[0]} = $t[0];
}
my %gen = ();
open my $pfh, '<', 'putclust.tsv' or die $!;
while(<$pfh>){
    chomp;
    my ($a, $k, $n) = split(/\t/, $_);
    print join("\t", $a, $k, $n, $taxof{$a})."\n";
    $gen{$taxof{$a}} += 1;
}
close $pfh;

foreach my $g (keys %gen){
    print STDERR join("\t", $gen{$g}, $g)."\n";
}
