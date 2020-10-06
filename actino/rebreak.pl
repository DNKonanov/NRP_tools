#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my %done = ();
foreach my $ens (glob("a*/ens.res.tsv")){
    open my $efh, '<', $ens or die $!;
    while(<$efh>){
	chomp;
	my($dom, $m, $spec) = split(/\t/, $_);
	$done{$dom} = 1;
    }
    close $efh;
}
print 83589-scalar(keys %done)." sequences remaining\n";
my $i = 1;
my $a = 9;
mkdir 'a'.$a unless(-d 'a'.$a);
open my $afh, '>', 'a'.$a."/adom.$a.faa" or die $!;
my $fa = new Bio::SeqIO(-file=>'adoms.faa', -format=>'fasta');
while(my $seq=$fa->next_seq){
    unless(exists $done{$seq->id}){
	if($i % 2950 == 0){
	    close $afh;
	    $a+=1;
	    mkdir 'a'.$a unless(-d 'a'.$a);
	    open $afh, '>', 'a'.$a."/adom.$a.faa" or die $!;
	}
	print $afh '>'.$seq->id."\n".$seq->seq."\n";
	$i+=1;
    }
}
close $afh;
