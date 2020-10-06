#!/bin/env perl

use strict;
use warnings;

my $tota = 20296;
chomp(my $n = `wc -l seqsim_matrix.txt`);
$n =~ s/^(\d+).+/$1/;
$n-=1;
my $nsum = ($n*($n+1))/2;
my $tsum = ($tota*($tota+1))/2;
my $togo = $tota-$n;
print sprintf("%.3f", $nsum/$tsum)." complete\t($nsum of $tsum)\t($togo remaining)\n";
