import argparse
import numpy as np
import scipy as sp

from detect_substrates import detect_substrates
from find_positions import gen_ex_adomains_file
from align import find_biosynthesis_clusters, custom_align

from subprocess import call
from Bio.SeqIO import parse
from Bio.Seq import translate, reverse_complement

parser = argparse.ArgumentParser()

parser.add_argument('-genome', type=str, required=True)
parser.add_argument('-SMILE', type=str, required=True)
parser.add_argument('-method', type=str, default='sandpuma')
parser.add_argument('-library', type=str, default='library.csv')

args = parser.parse_args()


print('Substrates detection... ')
substrates = detect_substrates(args.SMILE, args.library)

print()
print('Substrates: \n')
print(*substrates, sep=' -> ')
print()



call('./extract_adomains.sh {genome} '.format(genome=args.genome), shell=True)
print('Done!')

print('Find adomains positions...')
gen_ex_adomains_file('adom.faa', args.genome)


print('Substrates predictions...')
call('./predictnrps.sh {adomains} .'.format(adomains='adom.ex.faa'), shell=True)


clusters, cluster_positions  = find_biosynthesis_clusters('ind.res.tsv', eps=100000, min_samples=1)

print('Target cluster search...')
for i in clusters:
    pos = cluster_positions[i]
    print('position: ', pos)
    A = custom_align(substrates, clusters[i])
    print(A)
    print('\n')

print('Complete!')




