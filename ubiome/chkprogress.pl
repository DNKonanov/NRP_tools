#!/bin/env perl

use strict;
use warnings;

my $tot = 1955;
my @comp = split("\n", `wc -l */ens*v`);
my @log = split("\n", `wc -l */log`);
my $done = 0;
foreach(@comp){
	if($_ =~ /total/){
		$_ =~ s/\s*(\d+)\s*.+/$1/;
		$done = $_;
		last;
	}
}
my $te = 0;
foreach(@log){
	if($_ =~ /a\d+/){
		$_ =~ s/^\s*//;
		my ($e, $f) = split(/\s+/, $_);
		if($e > 0){
			print "$e ERRORS with $f\n";
			$te += 1;
		}
	}
}
print sprintf("%.2f", 100*$done/$tot) . "% completed";
print ' (no errors)' if($te == 0);
print "\n";
