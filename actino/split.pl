#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

my $fa = new Bio::SeqIO(-file=>'adoms.faa', -format=>'fasta');
my $i = 1;
my $a = 1;
open my $afh, '>', 'a'.$a."/adom.$a.faa" or die $!;
while(my $seq=$fa->next_seq){
    if($i % 10802 == 0){
	close $afh;
	$a+=1;
	open $afh, '>', 'a'.$a."/adom.$a.faa" or die $!;
    }
    print $afh '>'.$seq->id."\n".$seq->seq."\n";
    $i+=1;
}
close $afh;
