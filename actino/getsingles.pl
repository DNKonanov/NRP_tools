#!/bin/env perl

use strict;
use warnings;

my %seen = ();
open my $ffh, '<', 'f1_network.tsv' or die $!;
while(<$ffh>){
    chomp;
    my ($s, $t, $d) = split(/\t/, $_);
    $seen{$s} = 1;
    $seen{$t} = 1;
}
close $ffh;
my %bl = ();
open my $cfh, '<', 'cbreaks.list' or die;
while(<$cfh>){
    chomp;
    $bl{$_} = 1;
}
close $cfh;
open my $pfh, '<', 'putclust.tsv' or die $!;
my $singles = 0;
my $blrmsingles = 0;
while(<$pfh>){
    chomp;
    my ($c, $k, $num) = split(/\t/, $_);
    unless(exists $seen{$k}){
	$singles += 1;
	unless(exists $bl{$k}){
	    $blrmsingles += 1;
	}
    }
}
close $pfh;
print "Full singles:\t$singles\nBLRM singles:\t$blrmsingles\n";
