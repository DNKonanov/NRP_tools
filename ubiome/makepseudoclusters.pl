#!/bin/env perl

use strict;
use warnings;

## Read in data
my %e = ();
open my $efh, '<', 'all.ens.tsv' or die $!;
while(<$efh>){
    chomp;
    my($adom, $method, $spec) = split(/\t/, $_);
    my @a = split(/_/, $adom);
    my $source = join('_', @a[0..$#a-2]);
    $e{$source}{$a[-2]} = {
	'spec' => $spec,
	'end' => $a[-1]
    };
}
close $efh;

## Assign clusters
my %last = ();
my $spacingcutoff = 15000;
my $k = 0;
my %clust = ();
foreach my $source (sort keys %e){
    foreach my $start (sort {$a <=> $b} keys %{$e{$source}}){
	if(exists $last{$source}){
	    if($start - $last{$source} > $spacingcutoff){ ## outside cluster range
		$k += 1;
	    }
	}else{
	    $k += 1;
	}
	$e{$source}{$start}{'cluster'} = $k;
	$clust{$k}{'adoms'} += 1;
	$clust{$k}{'source'} = $source;
	$last{$source} = $e{$source}{$start}{'end'};
    }
}

## Dump results for each Adomain
open my $afh, '>', 'all.clust.tsv' or die $!;
foreach my $source (sort keys %e){
    foreach my $start (sort {$a <=> $b} keys %{$e{$source}}){
	print $afh join("\t", $source, $start, $e{$source}{$start}{'end'}, $e{$source}{$start}{'spec'}, 'k'.$e{$source}{$start}{'cluster'})."\n";
    }
}
close $afh;

## Dump results for each cluster
my $adomcutoff = 3;
open my $kfh, '>', 'clust.tsv' or die $!;
open my $cfh, '>', 'putclust.tsv' or die $!;
foreach my $knum (sort {$a <=> $b} keys %clust){
    print $kfh join("\t", $clust{$knum}{'source'}, 'k'.$knum, $clust{$knum}{'adoms'})."\n";
    print $cfh join("\t", $clust{$knum}{'source'}, 'k'.$knum, $clust{$knum}{'adoms'})."\n" if($clust{$knum}{'adoms'} >= $adomcutoff);
}
close $kfh;
close $cfh;
