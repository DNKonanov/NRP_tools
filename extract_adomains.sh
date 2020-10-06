#!/bin/sh

## usage:
##      ./extract_adomains.sh <nucl.fna>

## Grab the A-domains
echo "Identifying NRPS adenylation domains..."
perl bin/adomhmm.pl $1
