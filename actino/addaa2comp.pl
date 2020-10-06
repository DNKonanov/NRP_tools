#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my %seqof = ();
my $fa = new Bio::SeqIO(-file=>'corenrps.faa', -format=>'fasta');
while(my $seq = $fa->next_seq){
    $seqof{$seq->id} = $seq->seq;
}
while(<>){
    chomp;
    my ($comp, $bgc) = split(/\t/, $_);
    print join("\t", $comp, $bgc, $seqof{$bgc})."\n";
}
