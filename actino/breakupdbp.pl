#!/bin/env perl

use strict;
use warnings;

system("rm dbp/*");
open my $dfh, '<', 'all.adom.dbp' or die $!;
my $last = '';
my $ofh = undef;
while(<$dfh>){
    chomp;
    my @a = split(/\t/, $_);
    my @b = split(/_/, $a[0]);
    my $g = join('_', @b[0..$#b-3]);
    open $ofh, '>>', 'dbp/'.$g.'.dbp' or die $!;
    print $ofh "$_\n";
    close $ofh;
}
close $dfh;
