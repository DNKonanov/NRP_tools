#!/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;


#my $treef = new Bio::TreeIO(-file=>shift, -format=>'newick');
#my $outf = new Bio::TreeIO(-file=>'>multilocus_rn.ph', -format=>'newick');
#while(my $tree = $treef->next_tree){
#    foreach my $leaf ($tree->get_leaf_nodes){
#	my @a = split(/_/, $leaf->id);
#	$leaf->id($a[0]);
#    }
#    $outf->write_tree($tree);
#}
while(<>){
    chomp;
    $_ =~ s/([^_]+)_[^:]+(:)/$1$2/g;
    print "$_\n";
}
