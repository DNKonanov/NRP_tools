#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::TreeIO;

my $aamin = 360; ## typical nrpsA is 515AA, so this is about 30% trimming
foreach my $arg (@ARGV){
	my $stem = $arg;
	$stem =~ s/(.*)\.fa(a)?/$1/;
	system("mkdir $stem") unless(-d "$stem");
	## Initial alignment
	print "Aligning & Pruning";
	system("mafft --quiet --namelength 40 $arg > $stem.afa");
	my $totree = undef;
	my $lastafa = "$stem.afa";
	if( failgap("$stem.afa", $aamin) ){
		print '.';
		my $it = 1;
		my $tfi = $stem . '_trim' . "$it.faa";
		my $na = $tfi;
		$na =~ s/\.faa/\.afa/;
		## First trim
		trimmer("$stem.afa", $tfi);
		## Second alignment
		system("mafft --quiet --namelength 40 $tfi > $na");
		while( failgap($na, $aamin) ){
			print '.';
			my $la = $na;
			$it += 1;
			$tfi = $stem . '_trim' . "$it.faa";
			trimmer($la, $tfi);
			$na = $tfi;
			$na =~ s/\.faa/\.afa/;
			system("mafft --quiet --namelength 40 $tfi > $na");
		}
		$lastafa = $na;
		## Clustal format
		$na =~ s/\.afa/\.aln/;
		system("mafft --clustalout --quiet --namelength 40 $tfi > $na");
		$totree = $na;
	}else{
		## Clustal format
		system("mafft --clustalout --quiet --namelength 40 $arg > $stem.aln");
		$totree = "$stem.aln";
	}
	print "DONE!\n";
	system("clustalw -TREE -INFILE=$totree");
	## Clean up
	my $mvpref = $lastafa;
	$mvpref =~ s/\.afa//;
	my @mv = ('*.ph', $mvpref . '*');
	system("mv $_ $stem") foreach(@mv);
	my @rm = ($stem . '.a*', $stem . '_*');
	system("rm $_") foreach(@rm);
	my ($treef, @rest) = glob("$stem/*.ph");
	my $faa = $treef;
	$faa =~ s/\.ph$/\.faa/;
	my $fin = $faa;
	$fin =~ s/_trim\d+\.faa//;
	my %lastspec = ();
	my @f = ();
	## Read tree
	my $treein = new Bio::TreeIO(-file=>$treef, -format=>'newick');
	while(my $tree = $treein->next_tree){
		## Loop through leaves to get distances
		## Due to large processing time for big trees, walking out from each node
		## is done rather than generating all distances
		print "Generating distance metrics...";
		my %distfrom = ();
		my @leaves = $tree->get_leaf_nodes;
		for (my $o=0;$o<scalar(@leaves);$o+=1){
			my $oid = $leaves[$o]->id;
			my $osp = getspec($oid);
			my ($bdis, $fdis, $allowed) = (0,0,2); ## Disagreement counters, allow for 2 before break
			for(my $b=$o-1;$b>=0;$b=$b-1){ ## Walk backwards
				last unless(defined $leaves[$b]);
				my $bid = $leaves[$b]->id;
				next if(exists $distfrom{$oid}{$bid}); ## Already processed
				my $bsp = getspec($bid);
				my $d = $tree->distance(-nodes => [$leaves[$o], $leaves[$b]]);
				$distfrom{$oid}{$bid} = $d;
				$distfrom{$bid}{$oid} = $d;
				$bdis += 1 if($bsp ne $osp);
				last if($bdis >= $allowed);
			}
			for(my $f=$o+1;$f<scalar(@leaves);$f+=1){ ## Walk forwards
				last unless(defined $leaves[$f]);
				my $fid = $leaves[$f]->id;
				next if(exists $distfrom{$oid}{$fid}); ## Already processed
				my $fsp = getspec($fid);
				my $d = $tree->distance(-nodes => [$leaves[$o], $leaves[$f]]);
				$distfrom{$oid}{$fid} = $d;
				$distfrom{$fid}{$oid} = $d;
				$fdis += 1 if($fsp ne $osp);
				last if($fdis >= $allowed);
			}
		}
		print "DONE!\n";
		## Loop through leaves to assign groups
		print "Assigning groups...";
		my %group = ();
		my %added = ();
		my $g = 1;
		foreach my $leafout (@leaves){
			my $oid = $leafout->id;
			next if(exists $added{$oid}); ## Seen
			my @oidarr = split(/_+/, $oid);
			my $g2add = $oid;
			$added{$oid} = 1;
			my $num = 1;
			foreach my $iid (sort {$distfrom{$oid}{$a} <=> $distfrom{$oid}{$b}} keys %{$distfrom{$oid}}){
				my @iidarr = split(/_+/, $iid);
				last if($oidarr[-1] ne $iidarr[-1]);
				$g2add .= ',' . $iid;
				$added{$iid} = 1;
				$num += 1;
			}
			next unless($num>1);  ## Comment this out to include hits of only 1
			$group{$g}{'spec'} = $oidarr[-1];
			$group{$g}{'mem'} = $g2add;
			$g += 1;
		}
		print "DONE!\n";
		## Grab group members
		print "Retrieving sequences...";
		foreach my $r (sort {$a <=> $b} keys %group){
			my %get = ();
			my @m = split(/,/, $group{$r}{'mem'});
			foreach(@m){
				$get{$_} = 'not_def';
			}
			## Read faa
			my $fa = new Bio::SeqIO(-file=>$faa, -format=>'fasta');
			my @s = ();
			while(my $seq = $fa->next_seq){
				my $cleanid = $seq->id;
				$cleanid =~ s/[\(\)]/_/g;
				if(exists $get{$cleanid}){
					$get{$cleanid} = $seq->seq;
				}
			}
			## Print faa
			$lastspec{$group{$r}{'spec'}} += 1;
			my $ofn = $group{$r}{'spec'} . '_' . $lastspec{$group{$r}{'spec'}} . '.faa';
			$ofn =~ s/\|/-/g;
			open my $ofa, '>', $ofn or die $!;
			foreach my $s (keys %get){
				die "Issue with seq $s\n" if($get{$s} eq 'not_def');
				print $ofa '>' . "$s\n" . $get{$s} . "\n";
			}
			close $ofa;
			push @f, $ofn;
		}
		print "DONE!\n";
	}
	my $dir = $fin . '_groupfiles';
	system("rm -r $dir") if(-d $dir);
	system("mkdir $dir");
	print "Generating alignments and pHMMs...";
	foreach my $faafi (@f){
		my $pref = $faafi;
		$pref =~ s/\.faa//;
		chomp(my $seqs = `grep '>' $faafi|wc -l`);
		open my $clwo, '>', "$pref.clw" or die $!;
		## Generate alignment
		system("mafft --clustalout --quiet --namelength 40 $faafi > tmp.clw");
		## Format alignment
		open my $clwi, '<', 'tmp.clw' or die $!;
		while(<$clwi>){
			$_ = "CLUSTAL W multiple sequence alignment\n" if($_ =~ m/^CLUSTAL/);
			print $clwo $_;
		}
		close $clwi;
		system("rm tmp.clw");
		close $clwo;
		## Generate pHMM
		system("hmmbuild $pref.hmm $pref.clw > /dev/null");
		## Cleanup
		my @mv = ($faafi, "$pref.clw", "$pref.hmm");
		foreach my $m(@mv){
			system("mv $m $dir/");
		}
	}
	print "DONE!\n";
	## Create and press full pHMM
	my $db = $fin . '_nrpsA.hmmdb';
	print "Creating and pressing full pHMM into $db...";
	foreach(glob("$db*")){
		system("rm $_");
	}
	system("cat $dir/*.hmm > $db");
	system("hmmpress $db > /dev/null");
	print "DONE!\n";
}

sub getspec{
	my @idarr = split(/_+/, $_[0]);
	return $idarr[-1];
}

sub failgap{ ## gaps/ext test
	my ($aln, $min) = ($_[0], $_[1]);
	my ($tot, $ngaps, $cgaps) = (0,0,0);
	my $afa = new Bio::SeqIO(-file=>$aln, -format=>'fasta');
	while(my $seq = $afa->next_seq){
		$tot += 1;
		if($tot == 1){ #for the first seq, check against min
			my $s = $seq->seq;
			$s =~ s/-//g;
			return 0 if(length($s) <= $min);
		}
		if($seq->seq =~ m/^-/){
			$ngaps += 1;	
		}
		if($seq->seq =~ m/-$/){
			$cgaps += 1;
		}
	}
	if($ngaps == 0 && $cgaps == 0){
		return 0;
	}else{
		return 1;
	}
}

sub trimmer{
	my ($aln, $out) = @_;
	## Parse and trim
	my ($ftrim, $rtrim) = (0, 0);
	my %gaps = ();
	my ($last, $current, $alignment) = ('','','');
	my $afa = new Bio::SeqIO(-file=>$aln, -format=>'fasta');
	while(my $seq = $afa->next_seq){
		if($seq->seq =~ m/^-/){
			$ftrim = 1;
			$gaps{$seq->id}{'left'} = 1;
		}
		if($seq->seq =~ m/-$/){
			$rtrim = 1;
			$gaps{$seq->id}{'right'} = 1;
		}
	}
	if($ftrim == 0 && $rtrim == 0){ ## no trimming
		system("cp $aln $out");
		die "Died.  Issue: No trims found necessary in $aln.  Copied to $out\n";
	}else{
		$aln =~ s/\.afa/\.faa/;
		open my $tfh, '>', $out or die $!;      
		my $faa = new Bio::SeqIO(-file=>$aln, -format=>'fasta');
		if($ftrim == 0){ ## C trim only
			while(my $seq = $faa->next_seq){
				$seq = $seq->trunc(1, $seq->length - 1) unless(exists $gaps{$seq->id}{'right'});
				print $tfh '>' . $seq->id . "\n" . $seq->seq . "\n";
			}
		}elsif($rtrim == 0){ ## N trim only
			while(my $seq = $faa->next_seq){
				$seq = $seq->trunc(2, $seq->length) unless(exists $gaps{$seq->id}{'left'});
				print $tfh '>' . $seq->id . "\n" . $seq->seq . "\n";
			}
		}else{ ## Both trim
			while(my $seq = $faa->next_seq){
				$seq = $seq->trunc(1, $seq->length - 1) unless(exists $gaps{$seq->id}{'right'});
				$seq = $seq->trunc(2, $seq->length) unless(exists $gaps{$seq->id}{'left'});
				print $tfh '>' . $seq->id . "\n" . $seq->seq . "\n";
			}
		}
		close $tfh;
	}
}
