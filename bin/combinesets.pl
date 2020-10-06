#!/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

my $newadom = shift;

## Setup paths
my $basedir = abs_path($0);
$basedir =~ s/bin\/combinesets\.pl//;

## Get latest fullset
my $fs = undef;
foreach ( sort { $b cmp $a } glob("$basedir/flat/fullset*") ){
	$fs = $_;
	last;
}

if(-z $newadom){
	print STDERR "Nothing to combine.\n";
}else{
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	$year += 1900;
	$mon += 1;
	$mon = '0' . $mon if($mon < 10);
	$mday = '0' . $mday if ($mday < 10);
	my $cat = 'fullset' . $year . $mon . $mday . '.faa';
	system("cat $fs $newadom > $cat");
}
