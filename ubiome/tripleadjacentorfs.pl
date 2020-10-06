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
    if(scalar(keys %{$k2o{$k}}) == 3){
	my $specstr = '';
	my $orfstr = '';
	my @orfs = sort keys %{$k2o{$k}};
	my @o1 = split(/_/, $orfs[0]);
	my @o2 = split(/_/, $orfs[1]);
	my @o3 = split(/_/, $orfs[2]);
	if($o1[-1] eq $o2[-1] && $o1[-1] eq $o3[-1]){
	    if($o1[-1] eq 'fwd'){ ## forward
		if($o1[-2]+1 == $o2[-2] && $o2[-2]+1 == $o3[-2]){
		    $specstr = join('_', @{$orf2spec{$orfs[0]}}, @{$orf2spec{$orfs[1]}}, @{$orf2spec{$orfs[2]}});
		    $orfstr = join('-', @orfs);
		}
	    }else{ ## reverse
		if($o1[-2]-1 == $o2[-2] && $o2[-2]-1 == $o3[-2]){
		    $specstr = join('_', reverse(@{$orf2spec{$orfs[2]}}), reverse(@{$orf2spec{$orfs[1]}}), reverse(@{$orf2spec{$orfs[0]}}));
		    $orfstr = join('-', reverse(@orfs));
		}
	    }
	}
	$spec2orfs{$specstr}{$orfstr} = 1;
    }
}
foreach my $ss (keys %spec2orfs){
    my $sz = () = $ss =~ /_/g;
    print join("\t", $ss, $sz, scalar(keys %{$spec2orfs{$ss}}), join('|', sort keys %{$spec2orfs{$ss}}))."\n";
}
