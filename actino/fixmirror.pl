#!/bin/env perl

use strict;
use warnings;

my %all = ();
while(<>){
    next if($_ =~ m/^node1/);
    chomp;
    my ($n1, $n2, $dist) = split(/\t/, $_);
    next if($n1 eq $n2);
    if(exists $all{$n2}{$n1}){
	if($dist > $all{$n2}{$n1}){
	    $all{$n2}{$n1} = $dist;
	}
    }else{
	$all{$n1}{$n2} = $dist;
    }
}
foreach my $n1 (keys %all){
    foreach my $n2 (keys %{$all{$n1}}){
	print join("\t", $n1, $n2, $all{$n1}{$n2})."\n";
    }
}
