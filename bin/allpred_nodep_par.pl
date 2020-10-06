#!/bin/env perl

use strict;
use warnings;
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
my $dirpref = join('/', @dir[0..$#dir-1]);
chdir $dirpref;
$qfile = $dir[-1];
open my $tofh, '>', 'ind.res.tsv' or die $!;
open my $eofh, '>', 'ens.res.tsv' or die $!;
open my $pidfh, '>', 'pid.res.tsv' or die $!;
my %qfa = %{fasta2hash($qfile)};
foreach my $id (keys %qfa){
	## Make tmp seq query
	open my $tf, '>', "tmpt.faa" or die $!;
	print $tf '>'.$id.'_'.$wildcard."\n".$qfa{$id}."\n";
	close $tf;
	## Run prediCAT
	my ($pc, $pcf, $nndist, $nnscore, $snnscore) = treescanner($knownfaa, 'tmpt.faa', $wildcard);
	print $tofh join("\t", $id, 'prediCAT_MP', $pc) . "\n";
	$pcf = 'no_confident_result' if($snnscore < $snnthresh);
	print $tofh join("\t", $id, 'prediCAT_SNN', $pcf, $snnscore) . "\n";
	## Run ASM
	stachextract('tmpt.faa');
	my $asm = code2spec($knownasm, 'tmpt.stach.faa');
	print $tofh join("\t", $id, 'ASM', $asm) . "\n";
	## Run SVM
	my $svm = processsvm('tmpt.faa');
	print $tofh join("\t", $id, 'SVM', $svm) . "\n";
	## Run pHMM
	my $phmm = processphmm($knownfaa, 'tmpt.faa');
	print $tofh join("\t", $id, 'pHMM', $phmm) . "\n";
	## Run Ensemble
	my $pid = getpid('tmpt.faa');
	print $pidfh join("\t", $id, $pid) . "\n";
	my $ens = ensemble($pc, $pcf, $asm, $svm, $phmm, $max_depth, $min_leaf_support, $jackknife_data, $pid);
	print $eofh join("\t", $id, 'ENS', $ens) . "\n";
	## Cleanup
	#system("rm tmp*");
}
close $tofh;
close $eofh;
close $pidfh;
system("rm etmp*");

sub fasta2hash{
    my %s = ();
    my @id = ();
    my @seq = ();
    open my $qfh, '<', shift or die $!;
    while(<$qfh>){
	chomp;
	if($_ =~ m/^>(.+)/){
	    push @id, $1;
	}else{
	    my $lastid = scalar(@id)-1;
	    if(exists $seq[$lastid]){
		$seq[$lastid] .= $_;
	    }else{
		push @seq, $_;
	    }
	}
    }
    for(my $i=0;$i<scalar(@id);$i++){
	$s{$id[$i]} = $seq[$i];
    }
    return(\%s);
}

sub getfirstseq{
    my ($id, $seq) = (undef, undef);
    open my $qfh, '<', shift or die $!;
    while(<$qfh>){
	chomp;
	if($_ =~ m/^>(.+)/){
	    last if(defined($id));
	    $id = $1;
	}else{
	    $seq .= $_;
	}
    }
    close $qfh;
    return($id, $seq);
}

sub treescanner{
	my ($basefaa, $q, $wc) = @_;
	## Grab first seq
	my %bfa = %{fasta2hash($basefaa)};
	my ($fid, $fseq) = getfirstseq($basefaa);
	open my $ttf, '>', 'tmp.trim.faa' or die $!;
	print $ttf '>'.$fid."\n".$fseq."\n";
	close $ttf;
	system("cat $q >> tmp.trim.faa");
	## Align query to first seq
	system("mafft --quiet --namelength 60 --op 5 tmp.trim.faa > tmp.trim.afa");
	## Trim query
	my %pfa = %{fasta2hash('tmp.trim.afa')};
	my ($head, $tail) = (undef, undef);
	my ($id, $s) = (undef, undef);
	foreach my $qseq (keys %pfa){
	    if($qseq =~ m/$wc/){
		($id, $s) = ($qseq, $pfa{$qseq});
	    }else{
		($head, $tail) = ($pfa{$qseq}, $pfa{$qseq});
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
    my $tree = shift;
    my $py = $basedir.'/bin/treeparserscore.py '; ## THIS NEEDS TO BE MADE AVAILABLE TO RUNNING FROM DIFFERENT FOLDERS
    system("python2 $py $tree > tmp.tp.tsv");
    ## Parse calls
    open my $tp, '<', 'tmp.tp.tsv' or die $!;
    while(<$tp>){
	chomp;
	my ($i, $c, $atn, $fc, $nn, $nnsc, $snnsc) = split(/\t/, $_);
	if(defined $c && defined $fc){
	    close $tp;
	    return ($c, $atn, $fc, $nn, $nnsc, $snnsc);
	}
    }
    close $tp;
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
	my %fa = %{fasta2hash($in)};
	open my $ofh, '>', "$out.tmp" or die $!;
	## Loop through query seqs
	foreach my $seq (keys %fa){
	    ## Make alignment to example stachelhaus seeds
	    open my $qout, '>', 'squery.faa' or die $!;
	    print $qout '>'.$seq."\n".$fa{$seq}."\n";
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
	    my %nfa = %{fasta2hash('nquery.afa')};
	    my ($stp, $enp) = (undef,undef);
	    my $bflag = 0;
	    foreach my $nseq (keys %nfa){
		next unless($nseq eq 'phe_grsA');
		my ($st, $en) = (1, 106);
		my $w = 0;
		my @s = split('', $nfa{$nseq});
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
		$issue{$seq} = 2;
		next;
	    }
	    ## Loop through again and trim
	    open my $qtrim, '>', 'qtrim.afa' or die $!;
	    my @pos = ();
	    foreach my $nseq (keys %nfa){
		my $sub = substr($nfa{$nseq}, $stp-1, $enp-$stp+1);
		print $qtrim '>'.$nseq."\n$sub\n";
		## Extract stachelhaus positions from the alignment
		if($nseq eq 'phe_grsA'){
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
		    if($1 eq $seq){
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
			    $issue{$seq} = 1;
			}
		    }
		}
	    }
	    close $qtin;
	    #print STDERR $seq->id . " Stach extraction completed with no errors.\n";
	}
	close $ofh;
	## Final loop to remove issue seqs
	my %ofa = %{fasta2hash("$out.tmp")};
	open my $fin, '>', $out or die $!;
	foreach my $seq (keys %ofa){
	    print $fin '>'.$seq."\n".$ofa{$seq}."\n" unless(exists $issue{$seq});
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
	my %tfa = %{fasta2hash($trainstach)};
	foreach my $seq (keys %tfa){
	    next if($tfa{$seq} =~ m/-/);
	    my @id = split(/_+/, $seq);
	    $code2spec{$tfa{$seq}}{$id[-1]} += 1;
	}
	my %qfa = %{fasta2hash($querystach)};
	foreach my $seq (keys %qfa){
	    my $call = undef;
	    if(exists $code2spec{$qfa{$seq}}){
		my %l = ();
		foreach my $y (keys %{$code2spec{$qfa{$seq}}}){
		    my @w = split(/\|/, $y);
		    foreach my $r (@w){
			$l{$r} += 1;
		    }
		}
		$call = join('|', sort keys %l);
	    }else{
		my %score = ();
		my @q = split(//, $qfa{$seq});
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
    system("python2 $codepred $faa");
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
    chomp(my $ml = `python2 $basedir/bin/classifyp.py $qmat $md $msl`);
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
