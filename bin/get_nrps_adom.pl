#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

my ($tosearch, $setlist) = (shift,shift);

## Setup paths
my $basedir = abs_path($0);
$basedir =~ s/bin\/get_nrps_adom\.pl//;
my $q = $basedir . '/flat/query_adom.faa';

unless(-s "$tosearch.nrpsA.bp"){
	## Get A_Domains
	my $evalue = '1e-20';
	system("makeblastdb -in $tosearch -out $tosearch.db -dbtype prot 2>&1 > /dev/null");
	system("blastp -query $q -db $tosearch.db -outfmt 6 -out $tosearch.nrpsA.bp -max_target_seqs 10000000 -num_threads 8 -evalue $evalue");
}

open my $bp, '<', "$tosearch.nrpsA.bp" or die $!;
my %s2e = ();
while(<$bp>){
	chomp;
	my ($query, $hit, $pctid, $alen, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/, $_);
	$s2e{$hit}{$sstart}{'end'} = $send;
}
## Get mod numbers
foreach my $h ( keys %s2e){
	my $m = 1;
	foreach my $s (sort { $a <=> $b } keys %{$s2e{$h}} ){
		$s2e{$h}{$s}{'mod'} = $m;
		$m+=1;
	}
}
close $bp;
my %g2s = ();
open my $lfh, '<', $setlist or die $!;
while(<$lfh>){
	next if ($_ =~ m/^Cluster/);
	chomp;
	my ($bgc, $mod, $evid, $spec, $gene) = split(/\t/, $_);
	## Translate module into numeric
	if($mod eq 'x'){
		$mod = 1;
	}elsif($mod =~ m/\w?(\d+)/){
		$mod = $1;
	}
	my $head = join('_', $bgc, $gene);
	$g2s{$head}{$mod} = $spec;
}
## Take care of mod 0 issues by adding 1
## And not starting at 1 issues by beg at 1
foreach my $h (keys %g2s){
	if(exists $g2s{$h}{'0'}){
		my @ord = ();
		foreach my $old (sort { $a <=> $b} keys %{$g2s{$h}}){
			push @ord, $g2s{$h}{$old};
		}
		$g2s{$h} = undef;
		for(my $o=0;$o<scalar(@ord);$o+=1){
			my $n = $o+1;
			$g2s{$h}{$n} = $ord[$o];
		}
	}
	unless(exists $g2s{$h}{'1'}){
		my @ord = ();
		foreach my $old (sort { $a <=> $b} keys %{$g2s{$h}}){
			push @ord, $g2s{$h}{$old};
		}
		$g2s{$h} = undef;
		for(my $o=0;$o<scalar(@ord);$o+=1){
			my $n = $o+1;
			$g2s{$h}{$n} = $ord[$o];
		}
	}
}
my $outfi = $tosearch;
$outfi =~ s/\.faa/_adom\.faa/;
open my $aout, '>', $outfi or die $!;
my $nrpsfa = new Bio::SeqIO(-file=>$tosearch, -format=>'fasta');
while(my $seq = $nrpsfa->next_seq){
	foreach my $s (sort { $a <=> $b } keys %{$s2e{$seq->id}} ){
		my $a_dom = $seq->trunc($s, $s2e{$seq->id}{$s}{'end'});
		if(exists $g2s{$seq->id}{ $s2e{$seq->id}{$s}{'mod'} } ){
			print $aout '>' . join('_', $seq->id, $s, 'mod' . $s2e{$seq->id}{$s}{'mod'}, $g2s{$seq->id}{ $s2e{$seq->id}{$s}{'mod'} } ) . "\n" . $a_dom->seq . "\n";
		}else{
			print STDERR "Issue with mod $s2e{$seq->id}{$s}{'mod'} of " . $seq->id . "...skipping.\n";
		}
	}
}
close $aout;
system("rm $tosearch.db* 2> /dev/null");
