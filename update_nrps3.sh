#!/bin/sh

## usage:
##	./update_nrps3.sh <MIBiG_json_dir>

## Find all the post-khayatt MIBiG jsons with NRPSs
echo "Finding new NRPSs in MIBiG..."
python bin/mibig_id_new.py $1 > pk_nrps.txt
echo "Extracting information for validated modules..."
python bin/mibig_mktbl.py pk_nrps.txt > pk_set.tsv
echo "Getting full NRPS sequences..."
perl bin/get_nrps_prot.pl pk_set.tsv > pk_nrps.faa
echo "Extracting A-domains from gathered NRPSs..."
perl bin/get_nrps_adom.pl pk_nrps.faa pk_set.tsv
echo "Combining with most recent set..."
perl bin/combinesets.pl pk_nrps_adom.faa
echo "Cleaning specificities..."
perl bin/cleanspec.pl
echo "Generating pHMM db..."
perl bin/buildphmmdb.pl *_cl.faa > /dev/null
echo "Cleaning up intermediate files..."
rm pk_* 2> /dev/null
