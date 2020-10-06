#!/bin env perl

use strict;
use warnings;

foreach my $fna (glob("~/genomes/20160811*/fna/*.fna")){
    my @p = split(/\//, $fna);
    print "Processing $p[-1]...\n";
    system("./extract_adomains.sh $fna");
    if(-z 'adom.faa'){
	system("rm adom.faa");
    }else{
	system("mv adom.faa actino/$p[-1].nrps.faa");
    }
    print "DONE\n";
}
