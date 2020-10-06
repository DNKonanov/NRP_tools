#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

## Get pairwise PIDs for each A-domain
system("diamond makedb --in adoms.faa -d all.adom") unless(-e 'all.adom.dmnd');
system("diamond blastp -q adoms.faa -d all.adom -e 1 -f tab -o all.adom.dbp -p 10 -k 1000000") unless(-e 'all.adom.dbp');
#open my $dfh, '<', 'all.adom.dbp' or die $!;
#my %fid = ();
#while(<$dfh>){
#    chomp;
#    my ($query, $hit, $pctid, $alen, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/, $_);
#    $fid{$query}{$hit} = $pctid/100;
#}
#close $dfh;

## Read in the spec data
my %sp = ();
open my $sfh, '<', 'corenrps.specs.txt' or die $!;
while(<$sfh>){
    chomp;
    my($single, $verbose) = split(/\t/, $_);
    $sp{$verbose} = $single;
}
close $sfh;

## Make the clade config file
open my $ofh, '>', 'corenrps.clade.txt' or die $!;
my $index = 1;
print $ofh join("\t", 'index', 'clade_index', 'PhyloNode_index', 'function_annotation', 'clade_index_ori')."\n";
foreach my $spec (sort keys %sp){
    print $ofh join("\t", $index, $sp{$spec}, $index, $sp{$spec}, 'index'.$index)."\n";
    $index += 1;
}
close $ofh;

## Read in the orf data
my %kos = ();
my %nameof = ();
open my $pfh, '<', 'all.clust.tsv' or die $!;
while(<$pfh>){
    chomp;
    my($contig, $orf, $start, $end, $spec, $k) = split(/\t/, $_);
    push @{$kos{$k}{$orf}}, $sp{$spec};
    my $ind = scalar(@{$kos{$k}{$orf}}) - 1;
    $nameof{$orf}{$ind} = join('_', $orf, $start, $end);
}
close $pfh;

## Read in putclusts
my %input = ();
open my $ifh, '<', 'putclust.tsv' or die $!;
while(<$ifh>){
    chomp;
    my ($contig, $k, $num) = split(/\t/, $_);
    $input{$k} = 1;
}
close $ifh;

## Make the mock-fasta
open my $cfo, '>', 'corenrps.faa' or die $!;
my %adom = ();
foreach my $k (sort keys %kos){
    next unless(exists $input{$k});
    my $clusterstring = '';
    my $orfnum = 1;
    my $clusta = 1;
    foreach my $o (sort keys %{$kos{$k}}){
	my $orfa = 1;
	my @r = split(/_/, $o);
	my $short = join('_', @r[0..$#r-1]);
	if($o =~ m/fwd$/){
	    my $p = 0;
	    foreach my $dom (@{$kos{$k}{$o}}){
		#$clusterstring .= $orfnum.$dom; ## add orf number
		$clusterstring .= $dom; ## AA sequence only
		$adom{$nameof{$o}{$p}} = {
		    'orfdir' => $r[-1],
		    'orf' => $short,
		    'orfdom' => $orfa,
		    'clustdom' => $clusta,
		    'spec' => $dom,
		    'clust' => $k,
		    'orfnum' => $orfnum,
		};
		$orfa += 1;
		$clusta += 1;
		$p += 1;
	    }
	}else{
	    my $p = scalar(@{$kos{$k}{$o}}) - 1;
	    foreach my $dom (reverse @{$kos{$k}{$o}}){
		#$clusterstring .= $orfnum.$dom; ## add orf number
		$clusterstring .= $dom; ## AA sequence only
		$adom{$nameof{$o}{$p}} = {
		    'orfdir' => $r[-1],
		    'orf' => $short,
		    'orfdom' => $orfa,
		    'clustdom' => $clusta,
		    'spec' => $dom,
		    'clust' => $k,
		    'orfnum' => $orfnum
		};
		$orfa += 1;
		$clusta += 1;
		$p -= 1;
	    }
	}
	$orfnum += 1;
    }
    print $cfo '>'."$k\n$clusterstring\n";
}
close $cfo;



## Make the identity matrix
my $done = 0;
my @ord = sort keys %adom;
my $torun = 0;
for my $k (1..scalar(@ord)-1){
    $torun += scalar(@ord)-$k;
}
my %truedomof = ();
my %tde = ();
for my $r (@ord){
    my @id = split(/_/, $r);
    my $i = join('_', @id[0..$#id-2]);
    push @{$truedomof{$i}}, $id[-2];
    $tde{$i}{$id[-2]} = $id[-1];
}
my %orfdomof = ();
my %ode = ();
chomp(my $h = `LC_ALL=C grep '>' adoms.faa`);
my @a = split(/\n/, $h);
for my $d (@a){
    $d =~ s/^>//;
    my @id = split(/_/, $d);
    my $i = join('_', @id[0..$#id-2]);
    push @{$orfdomof{$i}}, $id[-2];
    $ode{$i}{$id[-2]} = $id[-1];
}

#print STDERR "$torun pid to lookup\n";
my $mfh = undef;
my %done = ();
if(-e 'seqsim_matrix.txt'){
    open my $ifh, '<', 'seqsim_matrix.txt' or die $!;
    while(<$ifh>){
	chomp;
	my ($k1, $k2, $c, $adom, $spec, $attr) = split(/\|/, $_);
	$done{$adom} = 1;
    }
    close $ifh;
    open $mfh, '>>', 'seqsim_matrix.txt' or die $!;
}else{
    open $mfh, '>', 'seqsim_matrix.txt' or die $!;
}
my $d = 0;
foreach my $row (@ord){
    chomp(my $date = `date`);
    $d += 1;
    print STDERR "$date\tProcessing $d of 18988\n";
    next if(exists $done{$row});
    print $mfh join('|',
	       $adom{$row}{'clust'},
	       $adom{$row}{'clust'},
	       $adom{$row}{'orf'},
	       $row,
	       $adom{$row}{'spec'},
	       'A'.$adom{$row}{'clustdom'}
	);
    my $matrix = '';
    my $kill = 0;
    foreach my $col (@ord){
	$kill = 1 if($row eq $col);
	if($kill == 0){
	    my $f = 0;
	    my @r = split(/_/, $row);
	    my $rorf = join('_', @r[0..$#r-2]);
	    my @c = split(/_/, $col);
	    my $corf = join('_', @c[0..$#c-2]);
	    my $rind = 0;
	    foreach (sort {$a <=> $b} @{$truedomof{$rorf}} ){
		last if($_ == $r[-2]);
		$rind += 1;
	    }
	    my $cind = 0;
	    foreach (sort {$a <=> $b} @{$truedomof{$corf}} ){
		last if($_ == $c[-2]);
		$cind += 1;
	    }
	    my @r2 = sort {$a <=> $b} @{$orfdomof{$rorf}};
	    my $rowrn = join('_', $rorf, $r2[$rind], $ode{$rorf}{$r2[$rind]});
	    my @c2 = sort {$a <=> $b} @{$orfdomof{$corf}};
	    my $colrn = join('_', $corf, $c2[$cind], $ode{$corf}{$c2[$cind]});
	    my @r3 = split(/_/, $rowrn);
	    my $rowgen = join('_', @r3[0..$#r3-3]);
	    chomp(my $grep = `LC_ALL=C grep $colrn dbp/$rowgen.dbp`);
	    my %fid = ();
	    if($grep ne ''){
		foreach my $ln (split(/\n/, $grep)){
		    my ($query, $hit, $pctid, $alen, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore) = split(/\t/, $ln);
		    $fid{$query}{$hit} = $pctid/100;
		}
	    }
	    if(exists $fid{$rowrn}{$colrn}){
		$f = $fid{$rowrn}{$colrn};
	    }#elsif(exists $fid{$colrn}{$rowrn}){
		#$f = $fid{$colrn}{$rowrn};
	    #}else{
		#die "COMPARISON NOT FOUND:\n$rowrn\n$colrn\n";
	    #}
	    $matrix .= ','.$f;
	    $done += 1;
	    #print STDERR '.' if($done % 10000 == 0);
	}else{
	    last;
	}
    }
    print $mfh "$matrix\n";
}
#print STDERR "\n";
close $mfh;
