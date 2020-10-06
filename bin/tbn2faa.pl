#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

## usage:
##	perl tbn2faa.pl <blast.tbn> <source.fna>

my ($tbn, $source) = (shift, shift);
my %s = ();
my $fa = new Bio::SeqIO(-file=>$source, -format=>'fasta');
while(my $seq = $fa->next_seq){
	$s{$seq->id} = $seq;
}
open my $tfh, '<', $tbn or die $!;
while(<$tfh>){
	chomp;
	my ($query, $hit, $pctid, $alen, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/, $_);
	if(exists $s{$hit}){
		my $ss = '';
		if($sstart > $send){
			my $swap = $sstart;
			$sstart = $send ;
			$send = $swap;
			$ss = $s{$hit}->trunc($sstart, $send)->revcom->translate;
		}else{
			$ss = $s{$hit}->trunc($sstart, $send)->translate;
		}
		print '>' . join('_', $hit, $sstart, $send) . "\n". $ss->seq . "\n";
	}else{
		die "$0 DIED.  Issue: $hit not found in $source\n";
	}
}
close $tfh;
