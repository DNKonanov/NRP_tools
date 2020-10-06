#!/bin/env perl

use strict;
use warnings;

my %ann = ();
open my $cfh, '<', 'all.clust.tsv' or die $!;
my $c = 0;
while(<$cfh>){
    chomp;
    my ($contig, $start, $end, $spec, $k) = split(/\t/, $_);
    $c += 1 if(exists $ann{$contig}{$c});
    $ann{$contig}{$c} = {
	'start'=> $start,
	'end'=> $end,
	'spec'=> $spec,
	'k'=> $k
    };
}
close $cfh;

foreach my $gff ( glob("prod/*.gff") ){
    my %orfat = ();
    open my $gfh, '<', $gff or die $!;
    while(<$gfh>){
	next if($_ =~ m/^#/);
	chomp;
	my ($contig, $source, $type, $start, $end, $score, $strand, $frame, $attr) = split(/\t/, $_);
	if(exists $ann{$contig}){
	    my $orf = $attr;
	    $orf =~ s/ID=\d+_(\d+).+/$1/;
	    $orf = join('_', $contig, $orf);
	    if($strand eq '+'){
		$orf .= '_fwd';
	    }else{
		$orf .= '_rev';
	    }
	    foreach my $p ($start..$end){
		$orfat{$contig}{$p} = $orf;
	    }
	}
    }
    close $gfh;
    foreach my $contig (sort keys %orfat){
	foreach my $dom (sort { $a <=> $b } keys %{$ann{$contig}} ){
	    my $orf = 'none';
	    foreach my $q ($ann{$contig}{$dom}{'start'}..$ann{$contig}{$dom}{'end'}){
		if(exists $orfat{$contig}{$q}){
		    $orf = $orfat{$contig}{$q};
		    last;
		}
	    }
	    $ann{$contig}{$dom}{'orf'} = $orf;
	}
    }
}

foreach my $contig (sort keys %ann){
    foreach my $dom (sort {$a <=> $b} keys %{$ann{$contig}}){
	unless(exists $ann{$contig}{$dom}{'orf'}){
	    $ann{$contig}{$dom}{'orf'} = 'none';
	}
	print join("\t", $contig, $ann{$contig}{$dom}{'start'}, $ann{$contig}{$dom}{'end'},
		   $ann{$contig}{$dom}{'spec'}, $ann{$contig}{$dom}{'k'}, $ann{$contig}{$dom}{'orf'})."\n";
    }
}
