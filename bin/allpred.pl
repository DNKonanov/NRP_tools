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
	my ($pc, $pcf, $nndist, $nnscore, $snnscore) = treescanner($knownfaa, 'tmpt.faa', $wildcard);
	print $tofh join("\t", $seq->id, 'prediCAT_MP', $pc) . "\n";
	$pcf = 'no_confident_result' if($snnscore < $snnthresh);
	print $tofh join("\t", $seq->id, 'prediCAT_SNN', $pcf, $snnscore) . "\n";
	## Run ASM
	stachextract('tmpt.faa');
	my $asm = code2spec($knownasm, 'tmpt.stach.faa');
	print $tofh join("\t", $seq->id, 'ASM', $asm) . "\n";
	## Run SVM
	my $svm = processsvm('tmpt.faa');
	print $tofh join("\t", $seq->id, 'SVM', $svm) . "\n";
	## Run pHMM
	my $phmm = processphmm($knownfaa, 'tmpt.faa');
	print $tofh join("\t", $seq->id, 'pHMM', $phmm) . "\n";
	## Run Ensemble
	my $pid = getpid('tmpt.faa');
	print $pidfh join("\t", $seq->id, $pid) . "\n";
	## NOTE: POTENTIALLY ADD IF STATEMENT HERE TO MAKE "NO_CALL" IF CALLS NOT MADE IN TOP 4
	##       DATA EXPLORATION NEEDED NEEDED FIRST, CURRENTLY NOT IMPLEMENTED
	my $ens = ensemble($pc, $pcf, $asm, $svm, $phmm, $max_depth, $min_leaf_support, $jackknife_data, $pid);
	print $eofh join("\t", $seq->id, 'ENS', $ens) . "\n";
	## Cleanup
	system("rm tmp*");
}
close $tofh;
close $eofh;
close $pidfh;
system("rm etmp*");

sub treescanner{
	my ($basefaa, $q, $wc) = @_;
	## Grab first seq
	my $bfa = new Bio::SeqIO(-file=>$basefaa, -format=>'fasta');
	my $seq = $bfa->next_seq;
	open my $ttf, '>', 'tmp.trim.faa' or die $!;
	print $ttf '>' . $seq->id . "\n" . $seq->seq . "\n";
	close $ttf;
	system("cat $q >> tmp.trim.faa");
	## Align query to first seq
	system("mafft --quiet --namelength 60 --op 5 tmp.trim.faa > tmp.trim.afa");
	## Trim query
	my $pfa = new Bio::SeqIO(-format=>'fasta', -file=>'tmp.trim.afa');
	my ($head, $tail) = (undef, undef);
	my ($id, $s) = (undef, undef);
	while(my $qseq = $pfa->next_seq){
		if($qseq->id =~ m/$wc/){
			($id, $s) = ($qseq->id, $qseq->seq);
		}else{
			($head, $tail) = ($qseq->seq, $qseq->seq);
			if($head =~ m/^-/){
				$head =~ s/^(-+).+/$1/;
			}else{
				$head = '';
			}
			if($tail =~ m/-$/){
				$tail =~ s/\w(-+)$/$1/;
			}else{
				$tail = '';
			}
		}
	}
	$head = length($head);
	$tail = length($tail);
	substr($s, -1, $tail) = '';
	substr($s, 0, $head) = '';
	$s =~ s/-//g;
	## Compile all seqs
	open my $newfaa, '>', 'tmp.tree.faa' or die $!;
	print $newfaa '>' . $id . "\n" . $s . "\n";
	close $newfaa;
	system("cat $basefaa >> tmp.tree.faa");
	system("perl -pi -e 's/[\(\)]/-/g' tmp.tree.faa");
	## Align all seqs, make tree
	system("mafft --quiet --namelength 60 tmp.tree.faa > tmp.tree.aln");
	system("FastTree -quiet < tmp.tree.aln > tmp.tree.ph 2>/dev/null");
	## Make calls
	chomp(my $cwd = `pwd -P`);
	my ($i, $c, $fc, $nn, $nnsc, $snnsc) = treeparserscore("$cwd/tmp.tree.ph");
	if(defined $c && defined $fc){
		return ($c, $fc, $nn, $nnsc, $snnsc);
	}
}

## prediCAT predictor
sub treeparserscore{
	my $zero_thresh = 0.005;
	## Set normalization twin pairs
	my @nppref = ('Q70AZ7_A3', 'Q8KLL5_A3');
	my @npnode = ();
	## Set normalized distance cutoff for nearest neighbor
	## (empirically derived default = 2.5)
	my $dcut = 2.5;
	## Read tree
	my $treef = $_[0];
	my $treein = new Bio::TreeIO(-file=>$treef, -format=>'newick');
	while(my $tree = $treein->next_tree){
	        my @query = ();
	        my %leaves = ();
	        ## Loop through leaves to find groups
	        my ($last, $group, $leafnum) = (undef, 1, 1);
	        for my $leaf ($tree->get_leaf_nodes){
	                if($leaf->id =~ m/^$nppref[0]/){
	                        $npnode[0] = $leaf;
	                }elsif($leaf->id =~ m/^$nppref[1]/){
	                        $npnode[1] = $leaf;
	                }
	                my @id = split(/_+/, $leaf->id);
	                if(defined $last){ ## Every pass except the first
	                        unless($last eq $id[-1] || $last eq $wildcard){ ## begin a new group
	                                $group += 1;
	                        }
	                }
	                $last = $id[-1];
	                $leaves{$leafnum} = {
	                        'group' => $group, 
	                        'id'    => $leaf->id,
	                        'spec'  => $id[-1],
	                        'node'  => $leaf
	                };
	                ## Record queries
	                push @query, $leafnum if($id[-1] eq $wildcard);
	                $leafnum += 1;
	        }
	        foreach my $q (@query){
	                ## Get distances to knowns
	                my %distfromq = ();
	                foreach my $n (keys %leaves){
	                        if($q != $n && $leaves{$n}{'spec'} ne $wildcard){
	                                $distfromq{$n} = $tree->distance(-nodes => [$leaves{$q}{'node'}, $leaves{$n}{'node'}]);
	                        }
	                }
	                ## Sort distances
	                my @o = sort {$distfromq{$a} <=> $distfromq{$b}} keys %distfromq;
	                ## Get zero distances
	                my @z = ();
	                #print STDERR "query $q\n";
	                foreach my $o (@o){
	                        #print STDERR "\t $distfromq{$o}\n";
	                        if($distfromq{$o} <= $zero_thresh){
	                                push @z, $o if($distfromq{$o} >= 0);
	                        }else{
	                                last;
	                        }
	                }
	                my $forcedpred = 'no_force_needed';
	                my $pred = 'no_pred_made';
	                if(scalar(@z) > 0){ ## query has zero length known neighbors
	                        $pred = $leaves{$z[0]}{'spec'}; 
	                }else{
	                        #check it out
	                        $pred = checkclade($q, $q-1, $q+1, $wildcard, $tree, %leaves);
	                        if($pred eq 'deeperdive'){
	                                ## deeper dive bested on closest 2 to deal with nodes at local extremes
	                                ## force a prediction if still none
	                                ($pred, $forcedpred) = deeperdive($q, $tree, $o[0], $o[1], %leaves);
	                        }
	                }
	                my $normdist = $tree->distance(-nodes => [$npnode[0], $npnode[1]]);
	                my $nn = sprintf("%.3f", $distfromq{$o[0]} / $normdist);
	                my ($nnscore, $snnscore) = (0,0);
	                if($nn < $dcut){
	                        $snnscore = getscore($dcut, $normdist, \%distfromq, \%leaves, @o);
	                        $nnscore = calcscore($dcut, $nn);
	                }
	                return ($leaves{$q}{'id'}, $pred, $forcedpred, $nn, $nnscore, $snnscore);
	        }
	}
}

## Check for bounded prediCAT specificities
sub checkclade{
	my ($query, $lo, $hi, $wc, $tree, %l) = @_;
	if(exists $l{$lo} && exists $l{$hi}){ ## Not first or last
		if($l{$lo}{'spec'} eq $wc){ ## lower bound is wildcard
			checkclade($query, $lo-1, $hi, $wc, %l);
		}elsif($l{$hi}{'spec'} eq $wc){ ## upper bound is wildcard
			checkclade($query, $lo, $hi+1, $wc, %l);
		}else{
			## Get the lca's descendants and check specs
			my $lca = $tree->get_lca(-nodes => [$l{$lo}{'node'},$l{$hi}{'node'}]);
			my ($spec, $pass) = (undef, 1);
			for my $child ($lca->get_all_Descendents){
				if($child->is_Leaf){
					my @id = split(/_/, $child->id);
					if(defined $spec){
						unless($id[-1] eq $spec || $id[-1] eq $wc){
							$pass = 0;
						}
					}else{
						$spec = $id[-1];
					}
				}
			}
			if($pass == 0){
				return 'deeperdive';
			}else{
				return $spec;
			}
		}
	}else{ ## First or last
		return 'deeperdive';
	}
}

## Deeper look at neighbors, if checkclade fails
sub deeperdive{
	my ($query, $tree, $n, $p, %l) = @_;
	## Want q to nearest dist to be less than nearest to 2nd-nearest dist
	my $q_to_n = $tree->distance(-nodes => [$l{$query}{'node'},$l{$n}{'node'}]);
	my $n_to_pn = $tree->distance(-nodes => [$l{$n}{'node'},$l{$p}{'node'}]);
	my $q_to_pn = $tree->distance(-nodes => [$l{$query}{'node'},$l{$p}{'node'}]);
	if($q_to_n < $n_to_pn && $l{$n}{'spec'} eq $l{$p}{'spec'}){
		return ($l{$n}{'spec'}, 'no_force_needed');
	}elsif($q_to_n == $q_to_pn && $l{$n}{'spec'} ne $l{$p}{'spec'}){
		return ('no_confident_result', 'no_confident_result');
	}else{
		my $parent = $l{$query}{'node'}->ancestor;
		my @sister = sort {$tree->distance(-nodes => [$l{$query}{'node'},$a]) <=> $tree->distance(-nodes => [$l{$query}{'node'},$b]) } $parent->get_all_Descendents;
		foreach my $sis (@sister){
			if($sis->is_Leaf && $sis->id ne $l{$query}{'id'}){
				my @fc = split(/_+/, $sis->id);
				return ('no_confident_result', $fc[-1]);
				last;
			}
		}
		return ('no_confident_result', 'no_confident_result');
	}
}

## Sum the SNN score
sub getscore{
	my ($scaleto, $nd, $dref, $lref, @o) = @_;
	my %dist2q = %$dref;
	my %leaf = %$lref;
	my $score = 0;
	my $nnspec = $leaf{$o[0]}{'spec'};
	foreach my $o (@o){
		my $curspec = $leaf{$o}{'spec'};
		if($nnspec eq $curspec){
			my $tmpscore = calcscore($scaleto, $dist2q{$o} / $nd );
			if($tmpscore > 0){
				$score += $tmpscore;
			}else{
				last;
			}
		}else{
			last;
		}
	}
	return sprintf("%.3f", $score);
}

## Individual SNN scores
sub calcscore{
	my ($scaleto, $distance) = @_;
	if($distance >= $scaleto){
		return 0;
	}else{
		return ($scaleto - $distance) / $scaleto;
	}
}

## Extract codes for the ASM
sub stachextract{
	my $in = $_[0];
	my $out = $in;
	$out =~ s/\.faa/\.stach\.faa/;
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

## Assign ASM code specificities
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

## Predict on SVM Models
sub processsvm{
	my $faa = $_[0];
	## tmp files to cleanup after processing
	my @tmp = (
		'nrpspredictor2_codes.txt',
		'input.sig',
		'muscle.fasta',
		'infile.fasta',
		'query.rep'
	);
	## Make sure models are current
	unless(-e "$basedir/dependencies/NRPSPredictor2/data/models/NRPS2_SINGLE_CLUSTER/[gly].mdl"){
		system("rm -r $basedir/dependencies/NRPSPredictor2/data/models/NRPS2_SINGLE_CLUSTER") if(-d "$basedir/dependencies/NRPSPredictor2/data/models/NRPS2_SINGLE_CLUSTER");
		system("cp -r $basedir/flat/svm/NRPS2_SINGLE_CLUSTER $basedir/dependencies/NRPSPredictor2/data/models/NRPS2_SINGLE_CLUSTER");
	}
	## Set up environment
	my $nrps2 = "$basedir/dependencies/NRPSPredictor2/NRPSpredictor2.sh";
	my $codepred = "$basedir/dependencies/NRPSPredictor2/nrpscodepred.py";
	## Get signatures
	system("python3 $codepred $faa");
	## Run SVMs
	system("$nrps2 $basedir -i input.sig -r query.rep -s 1 > /dev/null");
	## Get results
	open my $rep, '<', 'query.rep' or die $!;
	while(<$rep>){
		next if($_ =~ m/^#/);
		chomp;
		my @res = split(/\t/, $_);
		my ($seqid, $svmpred) = ($res[0], $res[6]);
		## Cleanup
		#foreach my $t (@tmp){
		#	system("rm $t") if(-e $t);
		#}
		return($svmpred);
	}
	close $rep;
}

sub ensemble{
	my ($pc, $pcf, $asm, $svm, $phmm, $md, $msl, $dfile, $pid) = @_;
	unless(-e 'etmp.specmap.tsv' && -e 'etmp.features.tsv' && -e 'etmp.labels.tsv'){
		## Read in the data
		my %da = ();
		my %allspec = ();
		open my $dfh, '<', $dfile or die $!;
		while(<$dfh>){
			chomp;
			next if ($_ =~ m/^shuffle/);
			my ($shuf, $jk, $query, $pid, $truespec, $called_spec, $method, $covered, $correct, $methshuf, $uname, $bin) = split(/\t+/, $_);
			unless(exists $da{$uname}){
				$da{$uname}{'true'} = $truespec;
				$da{$uname}{'pid'} = $pid;
				$da{$uname}{'shuf'} = $shuf;
				$da{$uname}{'jk'} = $jk;
				$da{$uname}{'query'} = $query;
				$da{$uname}{'bin'} = $bin;
			}
			$called_spec = 'no_call' if($covered eq 'N');
			$da{$uname}{'method'}{$method} = $called_spec;
			$allspec{$truespec} = -1;
			$allspec{$called_spec} = -1;
		}
		close $dfh;
		## Map specificities to integers
		my $id = 0;
		my %i2s = ();
		open my $sfh, '>', 'etmp.specmap.tsv' or die $!;
		foreach my $spec (sort keys %allspec){
			$allspec{$spec} = $id;
			print $sfh join("\t", $id, $spec) . "\n";
			$i2s{$id} = $spec;
			$id += 1;
		}
		close $sfh;
		## Create feature and label files
		my @m = ('prediCAT', 'forced_prediCAT_snn50', 'svm', 'stach', 'phmm');
		## FEATURE ORDER: pid', 'prediCAT', 'forced_prediCAT_snn50', 'svm', 'stach', 'phmm')
		open my $ffh, '>', 'etmp.features.tsv' or die $!;
		open my $lfh, '>', 'etmp.labels.tsv' or die $!;
		foreach my $uname (keys %da){
			foreach(@m){
				$da{$uname}{'method'}{$_} = 'no_call' unless(exists $da{$uname}{'method'}{$_});
			}
			print $lfh $allspec{$da{$uname}{'true'}} . "\n";
			print $ffh join("\t", $da{$uname}{'pid'},
				getmatrix($da{$uname}{'method'}{'prediCAT'}, "\t"),
				getmatrix($da{$uname}{'method'}{'forced_prediCAT_snn50'}, "\t"),
				getmatrix($da{$uname}{'method'}{'svm'}, "\t"),
				getmatrix($da{$uname}{'method'}{'stach'}, "\t"),
				getmatrix($da{$uname}{'method'}{'pHMM'}, "\t")
			) . "\n";
		}
		close $ffh;
		close $lfh;
	}
	my %i = ();
	open my $ifh, '<', 'etmp.specmap.tsv' or die $!;
	while(<$ifh>){
		chomp;
		my ($n, $sp) = split(/\t/, $_);
		$i{$n}=$sp;
	}
	close $ifh;
	my $qmat = join(' ', $pid,
		getmatrix($pc, ' '),
		getmatrix($pcf, ' '),
		getmatrix($svm, ' '),
		getmatrix($asm, ' '),
		getmatrix($phmm, ' ')
	);
	chomp(my $ml = `python $basedir/bin/classifyp.py $qmat $md $msl`);
	return $i{$ml};
}

sub getmatrix{
	my ($spec, $del) = @_;
	my $mat = '';
	my %i2s = ();
	open my $sfh, '<', 'etmp.specmap.tsv' or die $!;
	while(<$sfh>){
		chomp;
		my ($id, $sp) = split(/\t/, $_);
		$i2s{$id} = $sp;
	}
	close $sfh;
	for(my $i=0;$i<scalar(keys %i2s);$i+=1){
		if($mat eq ''){
			if($spec eq $i2s{$i}){
				$mat .= '1';
			}else{
				$mat .= '0';
			}
		}else{
			if($spec eq $i2s{$i}){
				$mat .= $del . '1';
			}else{
				$mat .= $del . '0';
			}
		}
	}
	return $mat;
}

sub getpid{
	my $query = shift;
	my $db = $knownfaa;
	$db =~ s/\.faa$/\.db/;
	unless(-e "$db.pin"){
		system("makeblastdb -in $knownfaa -out $db -dbtype prot 2>&1 > /dev/null");
	}
	system("blastp -query $query -db $db -outfmt 6 -out tmp.pid.bp -evalue 1e-10 -num_threads 10");
	my ($bit, $pid) = (0,0);
	open my $bfh, '<', 'tmp.pid.bp' or die $!;
	while(<$bfh>){
		chomp;
		my ($bquery, $hit, $pctid, $alen, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/, $_);
		if($bit < $bitscore){
			$pid = $pctid;
			$bit = $bitscore;
		}
	}
	close $bfh;
	system("rm tmp.pid.bp");
	return $pid;
}

sub processphmm{
	my ($trainfaa, $query) = @_;
	my @tp = split(/\//, $trainfaa);
	$tp[-1] =~ s/\.faa$//;
	my $hmmdb = join('/', @tp, $tp[-1].'_nrpsA.hmmdb');
	system("hmmscan -o tmp.hmmscan.out --tblout tmp.hmmtbl.out --noali $hmmdb $query");
	return "no_call" unless(-e 'tmp.hmmscan.out');
	open my $sfh, '<', "tmp.hmmtbl.out" or die "Died in hmmscanner: $!";
	while(<$sfh>){
		chomp;
		if($_ =~ m/^(\w+\S*)\s/){
			my @s = split(/_+/, $1);
			close $sfh;
			return join('_', @s[0..$#s-1]);
		}
	}
	close $sfh;
	return "no_call";
}
