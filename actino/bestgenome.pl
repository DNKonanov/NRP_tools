#!/bin/env perl

use strict;
use warnings;

my %best = ();
open my $nfh, '<', 'n50.tsv' or die $!;
while(<$nfh>){
    next if($_=~m/^#/);
    chomp;
    my ($n, @rest) = split(/\t/, $_);
    my @t = split(/_/, $rest[-1]);
    $t[0] =~ s/fna\///;
    if(exists $best{$t[0]}){
	if($best{$t[0]}{'contigs'} > $n){
	    $best{$t[0]}{'contigs'} = $n;
	    $best{$t[0]}{'file'} = $rest[-1];
	}
    }else{
	$best{$t[0]}{'contigs'} = $n;
	$best{$t[0]}{'file'} = $rest[-1];
    }
}
close $nfh;
open my $pfh, '<', 'genuscount.tsv' or die $!;
mkdir 'best' unless(-d 'best');
while(<$pfh>){
    chomp;
    my ($c, $g) = split(/\t/, $_);
    print join("\t", $g, $best{$g}{'contigs'}, $best{$g}{'file'})."\n";
    system("cp ~/genomes/20160811-publicactinos/$best{$g}{'file'} best/");
}
close $pfh;
