#!/bin/env perl

use strict;
use warnings;

## Read in data
print STDERR "Reading results...";
my %e = ();
my %topull = ();
open my $efh, '<', 'sandpuma.res.tsv' or die $!;
while(<$efh>){
    chomp;
    my($adom, $method, $spec) = split(/\t/, $_);
    my @a = split(/_/, $adom);
    my $source = join('_', @a[0..$#a-3]);
    my $orf = join('_', @a[0..$#a-2]);
    $e{$source}{$orf}{'spec'} = $spec;
    $topull{$orf} = 1;
}
close $efh;

## Read in prodigal headers
print STDERR "\nReading orf headers...";
my %coord = ();
foreach my $faa ( glob("~/genomes/20160811-publicactinos/prod/*faa") ){
    chomp(my $p = `grep '>' $faa`);
    my @header = split(/\n/, $p);
    foreach my $h (@header){
	$h =~ s/^>(.+)/$1/;
	my ($orf, $start, $end, $dir, $attr) = split(/\s#\s/, $h);
	next unless(exists $topull{$orf});
	if($dir == -1){
	    my $swap = $start;
	    $start = $end;
	    $end = $start;
	}
	$coord{$orf} = {
	    'start' => $start,
	    'end' => $end,
	};
    }
}

## Assign clusters
print STDERR "\nAssigning clusters...";
my %last = ();
my $spacingcutoff = 10000;
my $k = 0;
my %clust = ();
foreach my $source (sort keys %e){
    foreach my $cur (sort orfsort keys %{$e{$source}}){
	if(exists $last{$source}){
	    die $cur."\n" unless(exists $coord{$cur}{'start'});
	    if($coord{$cur}{'start'} - $last{$source} > $spacingcutoff){ ## outside cluster range
		$k += 1;
	    }
	}else{
	    $k += 1;
	}
	$e{$source}{$cur}{'cluster'} = $k;
	$clust{$k}{'adoms'} += 1;
	$clust{$k}{'source'} = $source;
	$last{$source} = $coord{$cur}{'end'};
    }
}

## Dump results for each Adomain
print STDERR "\nDumping results...";
open my $afh, '>', 'all.clust.tsv' or die $!;
foreach my $source (sort keys %e){
    foreach my $orf (sort orfsort keys %{$e{$source}}){
	print $afh join("\t", $source, $orf, $coord{$orf}{'start'}, $coord{$orf}{'end'}, $e{$source}{$orf}{'spec'}, 'k'.$e{$source}{$orf}{'cluster'})."\n";
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
print STDERR "\nDONE\n";

sub orfsort{
    my @aid = split(/_/, $a);
    my @bid = split(/_/, $b);
    return $aid[-1] <=> $bid[-1];
}
