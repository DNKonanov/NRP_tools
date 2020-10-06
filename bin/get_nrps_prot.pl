#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

## Setup paths
my $basedir = abs_path($0);
$basedir =~ s/bin\/get_nrps_prot\.pl//;
my $edirect = $basedir . '/dependencies/edirect';

## Get latest fullset
my $fs = undef;
foreach ( sort { $b cmp $a } glob("$basedir/flat/fullset*") ){
	$fs = $_;
	last;
}

## Read current list
my %curlist = ();
open my $cfh, '<', $fs or die $!;
while(<$cfh>){
	chomp;
	if($_ =~ m/^>(.+)/){
		my @a = split(/_+/, $1);
		$curlist{$a[0]} += 1 if($a[0] =~ m/^BGC/);
	}
}
close $cfh;

my %g2b = (); 
open my $lfh, '<', shift or die $!; 
while(<$lfh>){
	next if ($_ =~ m/^Cluster/);
	chomp;
	my ($bgc, $mod, $evid, $spec, $gene) = split(/\t/, $_);
	#unless($mod =~ m/\w?(\d+)/){
	#	print STDERR "Issue parsing $bgc:\n\tbgc:\t$bgc\n\tmod:\t$mod\n";
	#	next;
	#}
	next if(exists $curlist{$bgc});
	$g2b{$gene} = $bgc;
}
close $lfh;
foreach my $g (keys %g2b){
	fetch($g,1); 
	my $t = new Bio::SeqIO(-file=>'tmp.faa', -format=>'fasta');
	my $found = 0;
	while(my $seq = $t->next_seq){
		die "ERROR:\t$found sequences found for $g\n" if($found > 1);
		$found+=1;
		print '>' . join('_', $g2b{$g}, $g) . "\n" . $seq->seq . "\n";
	}
	system("rm tmp.err");
}
system("rm tmp.faa");

sub fetch{
	system("$edirect/esearch -db protein -query \"$_[0]\" | $edirect/edirect.pl -fetch -db protein -format fasta 2> tmp.err |grep . > tmp.faa");
	if(-s 'tmp.err'){
		if($_[1] <= 16){
			print STDERR "Retrying $_[0] ($_[1] second wait)\n";
			sleep $_[1];
			fetch($_[0], $_[1]*2);
		}else{
			print STDERR "Unable to retrieve $_[0] after 6 attempts. Consider adding manually.\n";
		}
	}
}
