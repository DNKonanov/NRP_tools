from Bio.SeqIO import parse
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--adoms', type=str, required=True)
parser.add_argument('--basedir', type=str, required=True)

args = parser.parse_args()


def ASM(faa):
    pass

basedir = args.basedir

def stachextract(infile):

    outfile = infile.replace('.faa', '.stach.faa')
    gapopen = 3.4
    seed = basedir + '/flat/seed.afa'

    properlen = 117

    grsAcode = {
        k: 1 for k in [4,5,8,47,68,70,91,99,100]
        }
    issue = {}
    fa = parse(infile, format='fasta')

    

    
    
