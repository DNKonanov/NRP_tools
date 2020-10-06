# **SANDPUMA**
*Specificity of AdenylatioN Domain Prediction Using Multiple Algorithms*

![Image Alt](./flat/sandpuma_logo.png)

Copyright © 2016 Marc Chevrette

If you find SANDPUMA useful in your research, please cite: Chevrette et al., 2017 https://doi.org/10.1093/bioinformatics/btx400

####################################################################

*This project is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.*

*This project is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.*

*You should have received a copy of the GNU General Public License
along with this program (filename LICENSE).  If not, see
<http://www.gnu.org/licenses/>.*

#####################################################################

# **This is the production repository of SANDPUMA and prediCAT.**

*The development repository can be found at bitbucket.org/chevrm/nrps2*

## **Code/method contributors:**
* Marc Chevrette (chevrm at gmail dot com):  Lead developer
* Fabian Aicheler:  SVM development
* Marnix Medema:  Developer/coordinator

## **Other key contributors:**
* Cameron Currie
* Oliver Kohlbacher

## **Prerequisite software and packages:**
* python (packages: json, glob, re, sys, os, csv, scipy, sklearn, numpy)
* perl (packages: Bio::SeqIO, Bio::TreeIO, Cwd 'abs_path')
* mafft
* FastTree
* ClustalW
* hmmscan (HMMER3)

#####################################################################

## **INSTALLATION:**
		## ensure all above dependencies are installed

	  	## Install dependencies listed in the apt repositories
	  	> sudo apt-get install python perl mafft ncbi-blast+ clustalw hmmer
		
		## Install python dependencies with pip. If pip not installed,
		## google how to set up pip
		> sudo pip install json glob re sys os csv scipy sklearn numpy

		## Install bioperl through the CPAN shell
		> sudo perl -MCPAN -e shell
		## Within the shell, enter below and choose defaults for all
		## questions
		   >> install Bio::SeqIO

		## Download the FastTree executable and add to your path
		> wget http://www.microbesonline.org/fasttree/FastTree
		> sudo chmod 777 FastTree
		> nano ~/.bashrc
		## Add FastTree to path
		## e.g.:
			export PATH=$PATH:/path/to/FastTree
		> source ~/.bashrc

		## set the SVM path

		## List the current path
		> pwd -P
		## Edit the SVM path
		> nano dependencies/NRPSPredictor2/NRPSpredictor2.sh
		## Change the path in variable NRPSBASEDIR to the
		## full path from pwd -P plus /dependencies/NRPSPredictor2
		## e.g.:
			export NRPS2BASEDIR=/home/mchevrette/git/sandpuma/dependencies/NRPSPredictor2

#####################################################################

## **Example Usages:**

### Update from a MIBiG json repository:
	./update_nrps3.sh <MIBiG_json_dir>

### Extract NRPS A-Domains from a nucleotide fasta:
	./extract_adomains.sh <nucl.fna>

### Run predictor on NRPS A-domains (protein):
	./predictnrps.sh <adomains.faa> <SANDPUMA base dir>