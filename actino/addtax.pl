#!/bin/env perl

use strict;
use warnings;

my %g = ();
my $tot = 0;
my %gof = ();
open my $pfh, '<', 'putclustgenus.tsv' or die $!;
while(<$pfh>){
    chomp;
    my ($contig, $k, $doms, $gen) = split(/\t/, $_);
    $gof{$k} = $gen;
    $g{$gen} += 1;
    $tot += 1;
}
close $pfh;
while(<>){
    if($_ =~ m/^Source/){
	print join("\t", 'Source', 'Target', 'dist', 'Source_Genus', 'Target_Genus')."\n";
    }else{
	chomp;
	my ($s, $t, $d) = split(/\t/, $_);
	print join("\t", $s, $t, $d, $gof{$s}, $gof{$t})."\n";
    }
}
foreach my $gen (sort {$g{$b} <=> $g{$a}} keys %g){
    print STDERR join("\t", $gen, $g{$gen}, sprintf("%.3f", $g{$gen}/$tot))."\n";
}
