#!/bin/env perl

use strict;
use warnings;


#my ($d1, $t1, $tg1) = processchk();
#sleep(300);
my ($d2, $t2, $tg2) = processchk();
#my $permin = ($d2-$d1)/5;
#my $mtg = $tg2/$permin;
print join("\t", 'ACTINO', sprintf("%.3f", $d2/$t2).' complete', "($d2 of $t2 completed, $tg2 remaining)")."\n";
#print 'Current rate: '.sprintf("%.1f", $permin).' per minute. '.sprintf("%.1f", $mtg/60).' hours remaining.'."\n";

sub processchk{
    chomp(my $tot = `grep '>' ~/git/sandpuma/actino/adoms.faa|wc -l`);
    chomp(my $d = `wc -l ~/git/sandpuma/actino/a*/ens.res.tsv`);
    my $done = 0;
    foreach my $l (split(/\n/, $d)){
	$l =~ s/^\s+//;
	my($c, $n) = split(/\s+/, $l);
	$done = $c if($n eq 'total');
    }
    my $togo = $tot-$done;
    return ($done, $tot, $togo);
}

