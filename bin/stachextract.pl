#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

## usage:       perl stachextract.pl fasta_of_adomains.faa

## Setup script path
my $seed_dir = abs_path($0);
$seed_dir =~ s/bin\/\S+\.pl/flat/;

my $in = shift;
my $out = $in;
$out =~ s/\.faa/\.stach\.faa/;
my $gapopen = 3.40;
my $seed = "$seed_dir/seed.afa";
my $properlen = 117;
my %grsAcode = (4=>1, 5=>1, 8=>1, 47=>1, 68=>1, 70=>1, 91=>1, 99=>1, 100=>1); ## Indexes for stachelhaus code in grsA WITHOUT gaps
my %issue = ();
my $fa = new Bio::SeqIO(-file=>$in, -format=>'fasta');
open my $ofh, '>', "$out.tmp" or die $!;
## Loop through query seqs
while(my $seq = $fa->next_seq){
	## Make alignment to example stachelhaus seeds
	open my $qout, '>', 'squery.faa' or die $!;
	print $qout '>' . $seq->id . "\n" . $seq->seq . "\n";
	close $qout;
	system("cat $seed >> squery.faa");
	system("mafft --quiet --op $gapopen --namelength 40 squery.faa > squery.afa");
	## Remove problematic newlines
	open my $sq, '<', 'squery.afa' or die $!;
	open my $nq, '>', 'nquery.afa' or die $!;
	my $current = '';
	while(<$sq>){
		chomp;
		if($_ =~ m/^>/){
			print $nq "\n" unless($current eq '');
			$current = $_;
			print $nq "$_\n";
		}else{
			print $nq $_;
		}
	}
	close $sq;
	close $nq;
	## Loop through alignment and get start and end sites based on guide seq
	my $nfa = new Bio::SeqIO(-file=>'nquery.afa', -format=>'fasta');
	my ($stp, $enp) = (undef,undef);
	my $bflag = 0;
	while(my $nseq = $nfa->next_seq){
		next unless($nseq->id eq 'phe_grsA');
		my ($st, $en) = (1, 106);
		my $w = 0;
		my @s = split('', $nseq->seq);
		foreach my $i (1..scalar(@s)){
			$w+=1 if($s[$i-1] =~ m/\w/);
			$stp = $i if($w == $st);
			$enp = $i if($w == $en);
		}
		unless(defined $stp && defined $enp){
			#print STDERR "ERROR: Issue with finding start/end seqs. Check MSA and consider higher gap open penalty\n";
			$bflag = 1;
		}
	}
	if($bflag == 1){
		$issue{$seq->id} = 2;
		next;
	}
	## Loop through again and trim
	$nfa = new Bio::SeqIO(-file=>'nquery.afa', -format=>'fasta');
	open my $qtrim, '>', 'qtrim.afa' or die $!;
	my @pos = ();
	while(my $nseq = $nfa->next_seq){
		my $sub = $nseq->subseq($stp, $enp);
		print $qtrim '>' . $nseq->id . "\n$sub\n";
		## Extract stachelhaus positions from the alignment
		if($nseq->id eq 'phe_grsA'){
			my $w = 0;
			my @s = split('', $sub);
			my %seen = ();
			foreach my $i (1..scalar(@s)){
				$w+=1 if($s[$i-1] =~ m/\w/);
				if(exists $grsAcode{$w}){
					push @pos, $i-1 unless(exists $seen{$w});
					$seen{$w} += 1;
				}
			}
		}
	}
	close $qtrim;
	## Parse and print out the code with error checking
	open my $qtin, '<', 'qtrim.afa' or die $!;
	my $pflag = 0;
	my %check = (
		'phe_grsA'	=> 'DAWTIAAIC',
		'asp_stfA-B2'	=> 'DLTKVGHIG',
		'orn_grsB3'	=> 'DVGEIGSID',
		'val_cssA9'	=> 'DAWMFAAVL'
	);
	my $cur = '';
	while (my $ln = <$qtin>){
		if($ln =~ m/^>(\S+)/){
			if($1 eq $seq->id){
				print $ofh $ln;
				$pflag = 1;
			}else{
				$cur = $1;
			}
		}else{
			chomp($ln);
			my $sc = '';
			foreach my $p (@pos){
				$sc .= substr($ln, $p, 1);
			}
			if($pflag == 1){
				print $ofh "$sc\n";
				$pflag = 0;
			}else{
				unless($check{$cur} eq $sc){
					#print STDERR "ERROR: Positive cntrl Stach codes do no match...\n" . $check{$cur} . " (correct $cur)\n$sc (returned $cur)\n";
					$issue{$seq->id} = 1;
				}
			}
		}
	}
	close $qtin;
	#print STDERR $seq->id . " Stach extraction completed with no errors.\n";
}
close $ofh;
## Final loop to remove issue seqs
my $ofa = new Bio::SeqIO(-file=>"$out.tmp", -format=>'fasta');
open my $fin, '>', $out or die $!;
while(my $seq = $ofa->next_seq){
	print $fin '>' . $seq->id . "\n" . $seq->seq . "\n" unless(exists $issue{$seq->id});
}
close $fin;
## Cleanup
my @toclean = ('nquery.afa', 'qtrim.afa', 'squery.afa', 'squery.faa', "$out.tmp");
system("rm $_") foreach(@toclean);
