import rdkit as rd
#import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Fingerprints.FingerprintMols import FingerprintMol
from rdkit.DataStructs import FingerprintSimilarity

names = [
    'FullName',
    'ShortName',
    'Letter',
    'MF',
    'SMILE',
]
library = pd.read_csv('library.csv', header=None, names=names, sep='\t')


def create_library(source):

    names = [
        'FullName',
        'ShortName',
        'Letter',
        'MF',
        'SMILE'
        ]
    file = pd.read_csv('library.csv', header=None, names=names, sep='\t')

    library = []
    for i in file.index:

        mol = Chem.MolFromSmiles(file.SMILE[i])
        library.append(
            (mol, file.ShortName[i].lower())
        )
    return library
    

def match_substrates(fragments, library):

    matches = []

    for f in fragments:
        for template in library:
            if (FingerprintSimilarity(
                    FingerprintMol(template[0]),
                    FingerprintMol(f)) > 0.9):
                matches.append(template[1])
    
    return matches
    

def hydrolise(mol):

    peptide_bond = Chem.MolFromSmiles('C(=O)NC')
    ester_bond = Chem.MolFromSmiles('C(=O)OC')
    peptide_ids = mol.GetSubstructMatches(peptide_bond)
    ester_ids = mol.GetSubstructMatches(ester_bond)

    nm = Chem.EditableMol(mol)
    
    bonds_ids = []
    for x, _, y, __ in peptide_ids:
        nm.RemoveBond(x,y)
        bonds_ids.append(nm.AddBond(x, nm.AddAtom(Chem.Atom('O')), Chem.BondType.SINGLE))

    for x, _, y, __ in ester_ids:
        nm.RemoveBond(x,y)
        bonds_ids.append(nm.AddBond(x, nm.AddAtom(Chem.Atom('O')), Chem.BondType.SINGLE))

    h_m = nm.GetMol()
    fragments = Chem.GetMolFrags(h_m, asMols=True)
    print()
    if len (fragments) == len(peptide_ids) + len(ester_ids):
        print('Cyclic structure!', end='')
    elif len (fragments) == len(peptide_ids) + len(ester_ids) - 1:
        print('Linear structure!')
    else:
        print('Unknown molecule topology!')
        
    return fragments


def detect_substrates(SMILE, library):
    substrates = []

    library = create_library(library)
    m = Chem.MolFromSmiles(SMILE)

    fragments = hydrolise(m)

    substrates = match_substrates(fragments, library)
    
    return substrates

