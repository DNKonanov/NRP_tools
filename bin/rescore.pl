#!/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';

## Get base dir
my @spath = split(/\//, abs_path($0));
my $basedir = join('/', @spath[0..$#spath-2]);

## Input files
my ($pidf, $indf, $ensf) = (shift, shift, shift);

## Define path cutoff
my $cutoff = 0.5;

## Read in decision tree node map
open my $nfh, '<', "$basedir/flat/nodemap.tsv" or die $!;
my %node = ();
while(<$nfh>){
    next if($_ =~ m/^#/);
    chomp;
    my($i, $parent, $parent_call, $dec, $thresh) = split(/\t/, $_);
    $node{$i} = {
	'parent' => $parent,
	'decision' => $dec,
	'parent_call' => $parent_call,
	'thresh' => $thresh
    };
}
close $nfh;

## Define paths
my @path = ();
foreach my $n (sort {$a <=> $b} keys %node){
    if($node{$n}{'decision'} =~ m/^LEAF_NODE/){
	my $p = $node{$n}{'parent'};
	my $traceback = $node{$p}{'decision'}.'%'.$node{$p}{'thresh'}.'-'.$node{$n}{'parent_call'}.'&LEAF_NODE-'.$n;
	while($p != 0){
	    my $t = '';
	    ($p, $t) = getparent($p);
	    $traceback = $t.'&'.$traceback;
	}
	push @path, $traceback;
    }
}

## Parse pid and individual results
my %p = ();
open my $pfh, '<', $pidf or die $!;
while(<$pfh>){
    chomp;
    my($query, $pid) = split(/\t/, $_);
    $p{$query}{'pid'} = $pid;
}
close $pfh;
open my $afh, '<', $indf or die $!;
my %snnof = ();
while(<$afh>){
    chomp;
    next if($_ =~ m/^shuffle/);
    my($query, $method, $called_spec, @rest) = split(/\t/, $_);
    $called_spec = cleannc($called_spec);
    $p{$query}{$method} = $called_spec;
    if($method eq 'prediCAT_SNN'){
	$snnof{$query} = $rest[0];
    }
}
close $afh;

## Loop through all results
my $totpassed = 0;
foreach my $q (keys %p){
    foreach my $pa (@path){
	my $pass = 1;
	my @dec = split(/&/, $pa);
	foreach my $d (@dec){
	    last if($d =~ m/LEAF_NODE/);
	    my($decision, $threshchoice) = split(/%/, $d);
	    my($thresh, $choice) = split(/-/, $threshchoice);
	    if($decision eq 'pid'){
		if($choice eq 'T'){ ## need greater than thresh to pass
		    if($thresh <= $p{$q}{'pid'}){
			$pass = 0;
			last;
		    }
		}else{ ## need less than or eq thresh to pass
		    if($thresh > $p{$q}{'pid'}){
			$pass = 0;
			last;
		    }
		}
	    }else{ ## not pid
		$decision = cleannc($decision);
		my @a = split(/_/, $decision);
		my $m = join('_', @a[0..$#a-1]);
		my $sp = $a[-1];
		if($choice eq 'T'){ ## less than 0.5, so NOT sp
		    die join("\t", $q, $m)."\n" unless(exists $p{$q}{$m});
		    if($p{$q}{$m} eq $sp){
			$pass = 0;
			last;
		    }
		}else{ ## matches sp
		    if($p{$q}{$m} ne $sp){
			$pass = 0;
			last;
		    }
		}
	    }
	}
	$p{$q}{'pass'} = $pass;
	$p{$q}{'path'} = $pa if($pass==1);
	$totpassed += $pass;
    }
}
#print STDERR "$totpassed passed total\n";

## Parse path accuracies
my %pathacc = ();
open my $ufh, '<', "$basedir/flat/traceback.tsv" or die $!;
while(<$ufh>){
    chomp;
    my ($pct, $n, $path) = split(/\t/, $_);
    $path =~ s/\S+&(LEAF_NODE-\d+)$/$1/;
    $pathacc{$path}{'pct'} = $pct;
    $pathacc{$path}{'n'} = $n;
}
close $ufh;

my %e = ();
my $changes = 0;
open my $efh, '<', $ensf or die $!;
while(<$efh>){
    chomp;
    my($query, $method, $called_spec)=split(/\t/, $_);
    my $path = $p{$query}{'path'};
    $path =~ s/\S+&(LEAF_NODE-\d+)$/$1/;
    my $acc = $pathacc{$path}{'pct'};
    my $n = $pathacc{$path}{'n'};
    if($acc < $cutoff){
	$called_spec = 'no_call';
	$changes += 1;
    }
    print join("\t", $query, 'SANDPUMA', $called_spec, 'SNN='.$snnof{$query}, 'PATHACC='."$acc")."\n";
}
close $efh;
print STDERR "$changes changes made\n";

sub getparent{
    my $n = shift;
    my $p = $node{$n}{'parent'};
    #my $c = $node{$n}{'parent_call'};
    my $c = $node{$p}{'thresh'}.'-'.$node{$n}{'parent_call'};# if($node{$n}{'decision'} eq 'pid');
    return ($p, $node{$p}{'decision'}.'%'.$c);
}
sub cleannc{
    my $s = shift;
    $s =~ s/_result//;
    $s =~ s/no_confident/nocall/;
    $s =~ s/N\/A/nocall/;
    $s =~ s/no_call/nocall/;
    return($s);
}
