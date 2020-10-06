#!/bin/env python

import json, glob, re, sys

bgclist = [line.rstrip('\n') for line in open( sys.argv[1] )]
print '\t'.join(['Cluster', 'Module', 'Evidence', 'Specificity', 'Gene'])
for bgcp in bgclist:
	with open( bgcp ) as data_file:
		data = json.load(data_file)
		if 'NRP' in data['general_params']:
			p = re.split('/', bgcp)
			p[-1] = p[-1].replace('.json', '')
			if 'nrps_genes' in data['general_params']['NRP']:
				for d in data['general_params']['NRP']['nrps_genes']:
					for m in d['nrps_module']:
						mod = str(m['module_nr'])
						a = m['a_substr_spec']
						spec, evid = 'none', 'none'
						if 'evidence_a_spec' in a:
							evid = a['evidence_a_spec']
						if 'nonprot_adom_spec' in a:
							spec = a['nonprot_adom_spec']
						elif 'prot_adom_spec' in a:
							spec = a['prot_adom_spec']
						if (evid != 'none') and (spec != 'none'):
							if evid != 'Sequence-based prediction':
								print '\t'.join([p[-1], str(mod), evid, spec, d['nrps_gene'] ])
							else:
								if 'compounds' in data['general_params']:
									print '\t'.join([p[-1], str(mod), evid, spec, d['nrps_gene'] ])
