#!/bin/env perl

use strict;
use warnings;

my %spec = ();
open my $sfh, '<', 'corenrps.specs.txt' or die $!;
while(<$sfh>){
    chomp;
    my ($short, $long) = split(/\t/, $_);
    $spec{$short} = $long;
}
close $sfh;
my %gof = ();
open my $gfh, '<', 'f1_tax_net.tsv' or die $!;
while(<$gfh>){
    next if ($_ =~ m/^Source/);
    chomp;
    my ($source, $target, $dist, $sg, $tg) = split(/\t/, $_);
    $gof{$source} = $sg;
    $gof{$target} = $tg;
}
close $gfh;
my ($lown, $medn, $hin) = (41, 43, 49);
my %comp = ();
open my $cfh, '<', 'f1comp_seq.tsv' or die $!;
while(<$cfh>){
    chomp;
    my ($c, $k, $s) = split(/\t/, $_);
    my $g = $gof{$k};
    $comp{$c}{'count'} += 1;
    $comp{$c}{'seq'}{$s} += 1;
    $comp{$c}{'gen'}{$g} += 1;
}
close $cfh;
foreach my $w (keys %comp){
    if($comp{$w}{'count'} == $lown){
	foreach my $s (keys %{$comp{$w}{'seq'}}){
	    my $full = '';
	    foreach my $sp (split(//, $s)){
		$full .= $spec{$sp} . ' ';
	    }
	    $full =~ s/\s$//;
	    print join("\t", $w, 'low', $full, $comp{$w}{'seq'}{$s})."\n";
	}
	foreach my $g (keys %{$comp{$w}{'gen'}}){
	    print STDERR join("\t", $w, 'low', $g, $comp{$w}{'gen'}{$g})."\n";
	}
    }elsif($comp{$w}{'count'} == $medn && $w eq 'comp18'){
	foreach my $s (keys %{$comp{$w}{'seq'}}){
	    my $full = '';
	    foreach my $sp (split(//, $s)){
		$full .= $spec{$sp} . ' ';
	    }
	    $full =~ s/\s$//;
	    print join("\t", $w, 'med', $full, $comp{$w}{'seq'}{$s})."\n";
	}
	foreach my $g (keys %{$comp{$w}{'gen'}}){
	    print STDERR join("\t", $w, 'med', $g, $comp{$w}{'gen'}{$g})."\n";
	}
    }elsif($comp{$w}{'count'} == $hin){
	foreach my $s (keys %{$comp{$w}{'seq'}}){
	    my $full = '';
	    foreach my $sp (split(//, $s)){
		$full .= $spec{$sp} . ' ';
	    }
	    $full =~ s/\s$//;
	    print join("\t", $w, 'hi', $full, $comp{$w}{'seq'}{$s})."\n";
	}
	foreach my $g (keys %{$comp{$w}{'gen'}}){
	    print STDERR join("\t", $w, 'hi', $g, $comp{$w}{'gen'}{$g})."\n";
	}
    }
}
