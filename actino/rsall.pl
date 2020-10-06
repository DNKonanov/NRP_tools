#!/bin env perl

use strict;
use warnings;

foreach my $i (1..20){
    print "a$i...\n";
    system("perl ../bin/rescore.pl a$i/pid.res.tsv a$i/ind.res.tsv a$i/ens.res.tsv > a$i.sp.tsv");
}
