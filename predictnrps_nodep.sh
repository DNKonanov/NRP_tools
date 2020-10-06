#!/bin/sh

## usage:
##      ./predictnrps.sh <adomains.faa> <SANDPUMA base dir>

## Run
perl $2/bin/allpred_nodep.pl $1
perl $2/bin/rescore.pl pid.res.tsv ind.res.tsv ens.res.tsv > sandpuma.tsv
