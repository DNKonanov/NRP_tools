#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $thresh = 10000;
my %sz = ();
foreach my $fna (glob("~/genomes/20160811*/fna/*.fna")){
    my $fa = new Bio::SeqIO(-file=>$fna, -format=>'fasta');
    while(my $seq = $fa->next_seq){
	$sz{$seq->id} = $seq->length;
    }
}
my %blacklist = ();
open my $afh, '<', 'all.clust.tsv' or die $!;
while(<$afh>){
    chomp;
    my ($contig, $orf, $start, $end, $spec, $clust) = split(/\t/, $_);
    if($start <= $thresh){
	$blacklist{$clust} += 1;
    }elsif($sz{$contig} - $end <= $thresh){
	$blacklist{$clust} += 1;
    }
}
close $afh;
foreach my $c (keys %blacklist){
    print "$c\n";
}
