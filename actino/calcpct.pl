#!/bin/env perl

use strict;
use warnings;

chomp(my $w = `wc -l seq*`);
foreach my $ln (split(/\n/, $w)){
    if($ln =~ m/(\d+)\stotal/){
	print sprintf("%.3f", $1/18988)." complete ($1 of 18988)\n";
    }
}
