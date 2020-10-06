#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $orig = shift;
my $pref = $orig;
$pref =~ s/\.faa$//;
my $break = 200;
my $i = 1;
my $b = 0;
my $fa = new Bio::SeqIO(-format=>'fasta', -file=>$orig);
my $bfh = undef;
while(my $seq=$fa->next_seq){
	if($i % $break == 0){
		close $bfh;
	}
	if($i % $break == 0 || $i == 1){
		$b+=1;
		open $bfh, '>', "$pref.$b.faa" or die $!;
	}
	print $bfh '>' . $seq->id . "\n" . $seq->seq . "\n";
	$i+=1;
}
