#!/bin/env perl

use strict;
use warnings;

my %strand = ();
foreach my $gff (glob("prod/*.rn.gff")){
    open my $gfh, '<', $gff or die $!;
    while(<$gfh>){
	next if($_ =~ m/^#/);
	chomp;
	my($contig, $source, $type, $start, $end, $score, $strand, $frame, $attr) = split(/\t/, $_);
	if($start > $end){
	    my $swap = $start;
	    $start = $end;
	    $end = $swap;
	}
	foreach my $p ($start..$end){
	    $strand{$source}{$p} = $strand;
	}
    }
    close $gfh;
}
