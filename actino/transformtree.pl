#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

## Read in the abbreviated AA list
my %s2v = ();
open my $mfh, '<', 'corenrps.specs.txt' or die $!;
while(<$mfh>){
    chomp;
    my ($s, $v) = split(/\t/, $_);
    $s2v{$s} = $v;
}
close $mfh;

## Read in the sequences and translate to verbose
my %k2v = ();
open my $cfh, '<', 'corenrps.tsv' or die $!;
while(<$cfh>){
    chomp;
    my ($k, $seq) = split(/\t/, $_);
    $seq =~ s/\|/\|\|\|/g;
    my $vstring = '';
    foreach my $c (split(//, $seq)){
	if($c eq '|'){
	    $vstring .= $c;
	}else{
	    unless($vstring eq '' || $vstring =~ m/\|$/){
		$vstring .= '-';
	    }
	    $vstring .= $s2v{$c};
	}
    }
    $k2v{$k} = $vstring;
}
close $cfh;

my $treei = new Bio::TreeIO(-file=>'mc_tree_upgmma_all_clade.nwk', -format=>'newick');
my $treeo = new Bio::TreeIO(-file=>'>mc_tree_trans.nwk', -format=>'newick');
while(my $tree = $treei->next_tree){
    foreach my $leaf ($tree->get_leaf_nodes){
	$leaf->id($leaf->id." ".$k2v{$leaf->id});
    }
    $treeo->write_tree($tree);
}
