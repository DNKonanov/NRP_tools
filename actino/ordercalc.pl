#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

## Get pairwise PIDs for each A-domain
#system("diamond makedb --in adoms.faa -d all.adom") unless(-e 'all.adom.dmnd');
#system("diamond blastp -q adoms.faa -d all.adom -e 1 -f tab -o all.adom.dbp -p 10 -k 1000000") unless(-e 'all.adom.dbp');
#open my $dfh, '<', 'all.adom.dbp' or die $!;
#my %fid = ();
#while(<$dfh>){
#    chomp;
#    my ($query, $hit, $pctid, $alen, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/, $_);
#    $fid{$query}{$hit} = $pctid/100;
#}
#close $dfh;

## Read in the spec data
my %sp = ();
open my $sfh, '<', 'corenrps.specs.txt' or die $!;
while(<$sfh>){
    chomp;
    my($single, $verbose) = split(/\t/, $_);
    $sp{$verbose} = $single;
}
close $sfh;

## Make the clade config file
open my $ofh, '>', 'corenrps.clade.txt' or die $!;
my $index = 1;
print $ofh join("\t", 'index', 'clade_index', 'PhyloNode_index', 'function_annotation', 'clade_index_ori')."\n";
foreach my $spec (sort keys %sp){
    print $ofh join("\t", $index, $sp{$spec}, $index, $sp{$spec}, 'index'.$index)."\n";
    $index += 1;
}
close $ofh;

## Read in the orf data
my %kos = ();
my %nameof = ();
open my $pfh, '<', 'all.clust.tsv' or die $!;
while(<$pfh>){
    chomp;
    my($contig, $orf, $start, $end, $spec, $k) = split(/\t/, $_);
    push @{$kos{$k}{$orf}}, $sp{$spec};
    my $ind = scalar(@{$kos{$k}{$orf}}) - 1;
    $nameof{$orf}{$ind} = join('_', $orf, $start, $end);
}
close $pfh;

## Read in putclusts
my %input = ();
open my $ifh, '<', 'putclust.tsv' or die $!;
while(<$ifh>){
    chomp;
    my ($contig, $k, $num) = split(/\t/, $_);
    $input{$k} = 1;
}
close $ifh;

## Make the mock-fasta
open my $cfo, '>', 'corenrps.faa' or die $!;
my %adom = ();
foreach my $k (sort keys %kos){
    next unless(exists $input{$k});
    my $clusterstring = '';
    my $orfnum = 1;
    my $clusta = 1;
    foreach my $o (sort keys %{$kos{$k}}){
	my $orfa = 1;
	my @r = split(/_/, $o);
	my $short = join('_', @r[0..$#r-1]);
	if($o =~ m/fwd$/){
	    my $p = 0;
	    foreach my $dom (@{$kos{$k}{$o}}){
		#$clusterstring .= $orfnum.$dom; ## add orf number
		$clusterstring .= $dom; ## AA sequence only
		$adom{$nameof{$o}{$p}} = {
		    'orfdir' => $r[-1],
		    'orf' => $short,
		    'orfdom' => $orfa,
		    'clustdom' => $clusta,
		    'spec' => $dom,
		    'clust' => $k,
		    'orfnum' => $orfnum,
		};
		$orfa += 1;
		$clusta += 1;
		$p += 1;
	    }
	}else{
	    my $p = scalar(@{$kos{$k}{$o}}) - 1;
	    foreach my $dom (reverse @{$kos{$k}{$o}}){
		#$clusterstring .= $orfnum.$dom; ## add orf number
		$clusterstring .= $dom; ## AA sequence only
		$adom{$nameof{$o}{$p}} = {
		    'orfdir' => $r[-1],
		    'orf' => $short,
		    'orfdom' => $orfa,
		    'clustdom' => $clusta,
		    'spec' => $dom,
		    'clust' => $k,
		    'orfnum' => $orfnum
		};
		$orfa += 1;
		$clusta += 1;
		$p -= 1;
	    }
	}
	$orfnum += 1;
    }
    print $cfo '>'."$k\n$clusterstring\n";
}
close $cfo;
my %a = ();
foreach my $sfile (glob("seq*txt")){
    open my $qfh, '<', $sfile or die $!;
    while(<$qfh>){
	chomp;
	my @b = split(/\|/, $_);
	$a{$b[3]} = $_;
    }
    close $qfh;
}

## Make the identity matrix

my @ord = sort keys %adom;
foreach my $row(@ord){
    print $a{$row}."\n";

}
