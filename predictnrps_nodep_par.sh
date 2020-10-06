#!/bin/sh

## usage:
##      ./predictnrps.sh <adomains.faa> <SANDPUMA base dir> <threads>

## Run
perl $2/bin/parbreakup.pl $1 $3
perl $2/bin/rescore.pl pid.res.tsv ind.res.tsv ens.res.tsv > sandpuma.tsv
