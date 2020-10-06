#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

my %spec = specload();
my %seq = ();
my $fn = undef;
foreach( sort { $b cmp $a } glob("fullset*.faa") ){
	next if($_ =~ m/_cl\.faa$/);
	$fn = $_;
	my $fa = new Bio::SeqIO(-file=>$fn, -format=>'fasta');
	while(my $seq=$fa->next_seq){
		my @header = split(/_+/, $seq->id);
		$header[-1] = lc($header[-1]);
		$header[-1] =~ s/^(\S+)/$1/;
		unless(exists $spec{$header[-1]}){ ## Not found in spec list
			system("clear");
			print STDERR "$header[-1] not found in list of specificities.\n";
			print STDERR 'What should the specificity of ' . $seq->id . " be added as? (how about $header[-1]?)\n";
			my @l = ();
			foreach my $k (keys %spec){
				if($spec{$k} eq 'keep'){
					push @l, $k;
				}else{
					push @l, $spec{$k};
				}
			}
			print STDERR "Current list:\n" . join("\t", sort @l) . "\n\n";
			chomp(my $new = <STDIN>);
			$new = $header[-1] if($new eq 'y');
			%spec = specload($new, $header[-1]);
			print STDERR "\n$header[-1] added to list of specificities as $new\n\n";
		}
		$header[-1] = $spec{$header[-1]} unless($spec{$header[-1]} eq 'keep');
		my $s = join('_', @header);
		$seq{$s} = $seq->seq;
	}
	last;
}
my $ofn = $fn;
$ofn =~ s/\.faa$/_cl\.faa/;
open my $ofh, '>', $ofn or die $!;
foreach(keys %seq){
	print $ofh  '>' . "$_\n" . $seq{$_} . "\n";
}
close $ofh;
system("rm $fn");

sub specload{ ## If(arg){adds to spec list}, loads spec list, & returns updated spec hash
	my $consol_dir = abs_path($0);
	$consol_dir =~ s/bin\/cleanspec\.pl/flat/;
	my $speclist = "$consol_dir/consol_speclist.txt";
	system("touch $speclist") unless(-e $speclist);
	if(defined $_[0] && defined $_[1]){
		open my $sfh, '>>', $speclist or die $!;
		if($_[0] eq $_[1]){
			print $sfh "$_[0]\n";
		}else{
			print $sfh join("\t", $_[1], $_[0]) . "\n";
		}
		close $sfh;
	}
	my %s = ();
	open my $si, '<', $speclist or die $!;
	while(<$si>){
		chomp;
		my ($sp, @ch) = split(/\t/, $_);
		if(scalar(@ch)>=1){
			foreach my $c (@ch){
				$s{$sp} = $c;
				$s{$c} = 'keep' unless(exists $s{$c});
			}
		}else{
			$s{$_} = 'keep';
		}
	}
	close $si;
	return %s;
}
