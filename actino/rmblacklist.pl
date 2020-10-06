#!/bin/env perl

use strict;
use warnings;

my %bl = ();
open my $cfh, '<', 'cbreaks.list' or die;
while(<$cfh>){
    chomp;
    $bl{$_} = 1;
}
close $cfh;
open my $nfh, '<', 'f1_tax_net.tsv' or die $!;
print join("\t", 'Source', 'Target', 'dist', 'Source_Genus', 'Target_Genus')."\n";
while(<$nfh>){
    chomp;
    my ($source, $target, @rest) = split(/\t/, $_);
    next if($source eq 'Source');
    next if(exists $bl{$source});
    next if(exists $bl{$target});
    print "$_\n";
}
close $nfh;
