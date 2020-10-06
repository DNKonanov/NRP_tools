#!/bin/env python

import json, glob, os, sys, re
from datetime import date

## Note, this will break on files of size 0 so runs
##	$ find mibigdir/ -size  0 -print0 |xargs -0 rm
## on the mibig json dir

os.system( 'find ' + sys.argv[1] + ' -size 0 -print0|xargs -0 rm 2> /dev/null' )
khayatt = date(2013,03,19) # set now to post-khayatt, can be changed in the future
dayssince = (date.today() - khayatt).days
for j in glob.glob( sys.argv[1] + '/*.json'):
	with open(j) as data_file:
		try:
			data = json.load(data_file)
			for c in data['general_params']['biosyn_class']:
				if (c == 'NRP') and ('publications' in data['general_params']):
					waszero = 0
					for pub in [data['general_params']['publications']]:
						for p in str(pub).split(','):
							if( not (re.search('\.', str(p) )) ): ## ignore non uid entries
								q = re.sub(r'\D', "", p)
								cmd = 'esearch -db pubmed -query "' + str(q) + '[uid]" | efilter -days ' + str(dayssince) + ' > tmp'
								os.system(cmd)
								with open('tmp') as tfh:
									content = tfh.readlines()
									count = re.search(r'<Count>(\d+)', str(content)).group(1)
									if(int(count) == 0):
										waszero = 1
					if(waszero == 0):
						print j
		except:
			sys.stderr.write( "Issue parsing " + j + "\n")
os.system('rm tmp')
