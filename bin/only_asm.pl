#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::TreeIO;
use Cwd 'abs_path';

## Get base dir
my @spath = split(/\//, abs_path($0));
my $basedir = join('/', @spath[0..$#spath-2]);

## Global Setting(s)
my $wildcard = 'UNK'; ## Default: 'UNK'

## PrediCAT Setting(s)
my $knownfaa = "$basedir/flat/fullset0_smiles.faa";
my $snnthresh = 0.5; ## Default: 0.5

## ASM Setting(s)
my $knownasm = "$basedir/flat/fullset0_smiles.stach.faa";

## Ensemble Setting(s)
my $max_depth = 40; ## Default 40
my $min_leaf_support = 10; ## Default 10
my $jackknife_data = "$basedir/flat/alldata.tsv";

## Run predictors
my $qfile = shift;
my @dir = split(/\//, $qfile);
my $fa = new Bio::SeqIO(-file=> $qfile, -format=> 'fasta');
open my $tofh, '>', 'ind.res.tsv' or die $!;
open my $eofh, '>', 'ens.res.tsv' or die $!;
open my $pidfh, '>', 'pid.res.tsv' or die $!;


while(my $seq = $fa->next_seq){
	## Make tmp seq query
	open my $tf, '>', "tmpt.faa" or die $!;
	print $tf '>' . $seq->id . '_' . $wildcard . "\n" . $seq->seq . "\n";
	close $tf;
	## Run prediCAT
	## my ($pc, $pcf, $nndist, $nnscore, $snnscore) = treescanner($knownfaa, 'tmpt.faa', $wildcard);
	## print $tofh join("\t", $seq->id, 'prediCAT_MP', $pc) . "\n";
	## $pcf = 'no_confident_result' if($snnscore < $snnthresh);
	## print $tofh join("\t", $seq->id, 'prediCAT_SNN', $pcf, $snnscore) . "\n";
	## Run ASM
	stachextract('tmpt.faa');
	my $asm = code2spec($knownasm, 'tmpt.stach.faa');
	print $tofh join("\t", $seq->id, 'ASM', $asm) . "\n";
	## Run SVM
	## my $svm = processsvm('tmpt.faa');
	## print $tofh join("\t", $seq->id, 'SVM', $svm) . "\n";
	## Run pHMM
	## my $phmm = processphmm($knownfaa, 'tmpt.faa');
	## print $tofh join("\t", $seq->id, 'pHMM', $phmm) . "\n";
	## Run Ensemble
	## my $pid = getpid('tmpt.faa');
	## print $pidfh join("\t", $seq->id, $pid) . "\n";
	## NOTE: POTENTIALLY ADD IF STATEMENT HERE TO MAKE "NO_CALL" IF CALLS NOT MADE IN TOP 4
	##       DATA EXPLORATION NEEDED NEEDED FIRST, CURRENTLY NOT IMPLEMENTED
	## my $ens = ensemble($pc, $pcf, $asm, $svm, $phmm, $max_depth, $min_leaf_support, $jackknife_data, $pid);
	## print $eofh join("\t", $seq->id, 'ENS', $ens) . "\n";
	## Cleanup
	system("rm tmp*");
}


sub code2spec{
	my ($trainstach, $querystach) = ($_[0], $_[1]);
	return('no_call') if(-z $querystach);
	my %code2spec = ();
	my $tfa = new Bio::SeqIO(-file=>$trainstach, -format=>'fasta');
	while(my $seq = $tfa->next_seq){
		next if($seq->seq =~ m/-/);
		my @id = split(/_+/, $seq->id);
		$code2spec{$seq->seq}{$id[-1]} += 1;
	}
	my $qfa = new Bio::SeqIO(-file=>$querystach, -format=>'fasta');
	while(my $seq = $qfa->next_seq){
		my $call = undef;
		if(exists $code2spec{$seq->seq}){
			my %l = ();
			foreach my $y (keys %{$code2spec{$seq->seq}}){
				my @w = split(/\|/, $y);
				foreach my $r (@w){
					$l{$r} += 1;
				}
			}
			$call = join('|', sort keys %l);
		}else{
			my %score = ();
			my @q = split(//, $seq->seq);
			foreach my $k (keys %code2spec){
				my @c = split(//, $k);
				my $cor = 0;
				for(my $i=0;$i<scalar(@q);$i+=1){
					$cor += 1 if($q[$i] eq $c[$i]);
				}
				push @{$score{$cor}}, $k;
			}
			foreach my $s (sort { $b <=> $a } keys %score){
				if($s >= 7){ ## at least 7 of 9 correct
					my %l = ();
					foreach my $k (@{$score{$s}}){
						foreach my $y (keys %{$code2spec{$k}}){
							my @w = split(/\|/, $y);
							foreach my $r (@w){
								$l{$r} += 1;
							}
						}
					}
					$call = join('|', sort keys %l);
				}else{
					$call = 'no_call';
				}
				last;
			}
		}
		return($call);
	}
}

sub stachextract{
	my $in = $_[0];
	print($in);
	
	my $out = $in;
	$out =~ s/\.faa/\.stach\.faa/;
	print($out);
	print('\n')
	my $gapopen = 3.40;
	my $seed = "$basedir/flat/seed.afa";
	die "ERROR: Unable to locate seed.afa for ASM alignment.\n" unless(-e $seed);
	my $properlen = 117;
	## Indicies for stachelhaus code in grsA WITHOUT gaps
	my %grsAcode = (4=>1, 5=>1, 8=>1, 47=>1, 68=>1, 70=>1, 91=>1, 99=>1, 100=>1);
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
			'phe_grsA'      => 'DAWTIAAIC',
			'asp_stfA-B2'   => 'DLTKVGHIG',
			'orn_grsB3'     => 'DVGEIGSID',
			'val_cssA9'     => 'DAWMFAAVL'
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
}
