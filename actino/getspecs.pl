#!/bin/env perl

use strict;
use warnings;

#my %sp = ();
#open my $sfh, '<', 'corenrps.specs.txt' or die $!;
#while(<$sfh>){
#    chomp;
#    my($single, $verbose) = split(/\t/, $_);
#    $sp{$verbose} = $single;
#}
#close $sfh;

my %s = ();
#my %kos = ();
open my $pfh, '<', 'all.clust.tsv' or die $!;
while(<$pfh>){
    chomp;
    my($contig, $orf, $start, $end, $spec, $k) = split(/\t/, $_);
    #push @{$kos{$k}{$orf}}, $sp{$spec};
    $s{$spec} += 1;
}
close $pfh;

#my %clust = ();
#foreach my $k (sort keys %kos){
#    my $clusterstring = '';
#    foreach my $o (sort keys %{$kos{$k}}){
#	if($o =~ m/fwd$/){
#	    $clusterstring .= '|'.join('', @{$kos{$k}{$o}});
#	}else{
#	    $clusterstring .= '|'.join('', reverse @{$kos{$k}{$o}});
#	}
#   }
#    $clusterstring =~ s/^\|//;
#    my @byorf = split(/\|/, $clusterstring);
#    $clusterstring = join('|', sort @byorf);
#    $clust{$clusterstring} += 1;
#    print join("\t", $k, $clusterstring)."\n";
#}

foreach my $s (sort keys %s){
    print STDERR "$s\n";
}

#foreach my $c (sort {$clust{$b} <=> $clust{$a}} keys %clust){
#    print STDERR join("\t", $clust{$c}, $c)."\n";
#}
