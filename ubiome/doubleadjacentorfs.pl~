#!/bin/env perl

use strict;
use warnings;

my %k2o = ();
my %orf2spec = ();
open my $pfh, '<', 'putorf.tsv' or die $!;
chomp(my @put = <$pfh>);
close $pfh;
foreach my $ln (@put){
    my ($contig, $start, $end, $spec, $k, $orf) = split(/\t/, $ln);
    $k2o{$k}{$orf} += 1;
    push @{$orf2spec{$orf}}, $spec;
}
my %spec2orfs = ();
foreach my $ln (@put){
    my ($contig, $start, $end, $spec, $k, $orf) = split(/\t/, $ln);
    if(scalar(keys %{$k2o{$k}}) == 1){
	my @o = split(/_/, $orf);
	my $specstr = '';
	if($o[-1] eq 'fwd'){
	    $specstr = join('_', @{$orf2spec{$orf}});
	}else{
	    $specstr = join('_', reverse(@{$orf2spec{$orf}}));
	}
	$spec2orfs{$specstr}{$orf} = 1;
    }
}
foreach my $ss (keys %spec2orfs){
    my $sz = () = $ss =~ /_/g;
    print join("\t", $ss, $sz, scalar(keys %{$spec2orfs{$ss}}), join('|', sort keys %{$spec2orfs{$ss}}))."\n";
}
