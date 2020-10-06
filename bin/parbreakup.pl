#!/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';

## Get base dir
my @spath = split(/\//, abs_path($0));
my $basedir = join('/', @spath[0..$#spath-2]);

my ($basefaa, $subfiles) = (shift, shift);
my @d = split(/\//, $basefaa);
my $pref = $d[-1];
$pref =~ s/\.faa$//;
my %bf = %{fasta2hash($basefaa)};
my $numseqs = scalar(keys %bf);
my $seqperf = int($numseqs / $subfiles)+1;
my $i = 0;
my $f = 1;
my $ffh = undef;
my @subf = ();
foreach my $seq (keys %bf){
    if($i == 0){
	open $ffh, '>', "$pref.$f.faa" or die $!;
	push @subf, "$pref.$f.faa";
	print $ffh '>'.$seq."\n".$bf{$seq}."\n";
	$i += 1;
    }
    elsif($i < $seqperf){
	print $ffh '>'.$seq."\n".$bf{$seq}."\n";
	$i += 1;
    }else{
	close $ffh;
	$f += 1;
	open $ffh, '>', "$pref.$f.faa" or die $!;
	push @subf, "$pref.$f.faa";
	print $ffh '>'.$seq."\n".$bf{$seq}."\n";
	$i = 1;
    }
}
close $ffh;
my $a = join(' ', @subf);
system("python2 $basedir/bin/par.py $a");
my @tocat = ('ind.res.tsv', 'ens.res.tsv', 'pid.res.tsv', 'query.rep');
foreach my $t (@tocat){
    system("rm $t") if(-e $t);
}
foreach my $s (@subf){
    my $p = $s;
    $p =~ s/\.faa$//;
    foreach my $t (@tocat){
	system("cat $p/$t >> $t");
    }
    ## Cleanup
    system("rm -r $p $s");
}

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
