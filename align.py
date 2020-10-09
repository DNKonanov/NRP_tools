import pandas as pd
from sklearn.cluster import DBSCAN
import string
import os
from skbio.util import classproperty
from skbio.sequence import GrammaredSequence
from skbio.alignment import global_pairwise_align, make_identity_substitution_matrix



letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

def find_biosynthesis_clusters(
        predicted_adoms,
        eps=100000,
        min_samples=1):
    predictions = pd.read_csv(
        predicted_adoms,
        sep='\t',
        header=None,
        names=[
            'header',
            'method',
            'result'
            ]
        )

    ASM = predictions.loc[predictions.method=='ASM']
    adoms = {}

    for i in ASM.index:
        h = predictions.header[i]
        coord = h.split('*')[-2].split(':')
        start, end = int(coord[0]), int(coord[1])
        strain = h.split('*')[-1]

        if (start, end) not in adoms:
            adoms[(start, end)] = []
        adoms[(start, end)].append(predictions.result[i])

    locs = list(adoms.keys())
    locs.sort()

    db_model = DBSCAN(eps=eps, min_samples=min_samples)

    X = []
    predicted_adoms = []
    for l in locs:
        X.append([l[0]])
        predicted_adoms.append(adoms[l][0])

    db_model.fit(X)
    clusters = {}
    cluster_positions = {}
    labels = db_model.labels_

    for i in range(len(labels)):
        if labels[i] not in clusters:
            clusters[labels[i]] = []
            cluster_positions[labels[i]] = []
        clusters[labels[i]].append(predicted_adoms[i])
        cluster_positions[labels[i]].append(locs[i])
        
    return clusters, cluster_positions
    

def generate_hash(target, cluster):

    target_set = set(target)
    cluster_set = set(cluster)

    intersect = target_set.intersection(cluster_set)

    hashes = {}
    letters_cnt = 0

    for i in intersect:
        hashes[i] = letters[letters_cnt]
        letters_cnt += 1

    for i in target_set:
        if i in hashes:
            continue
        hashes[i] = letters[letters_cnt]
        letters_cnt += 1

    for i in cluster_set:
        if i in hashes:
            continue
        hashes[i] = letters[letters_cnt]
        letters_cnt += 1
    return hashes

    

def code_seq(seq, hashes):
    new_seq = ''
    for i in seq:
        new_seq += hashes[i]
    return new_seq


def generate_combinations(seq):
    combinations = [seq]
    return combinations

def custom_align(target, cluster):
    hashes = generate_hash(target, cluster)
    coded_target = code_seq(target, hashes)
    coded_cluster = code_seq(cluster, hashes)

    print('Coded cluster: ', coded_cluster)

    class CustomSequence(GrammaredSequence):
        @classproperty
        def degenerate_map(cls):
            return {}

        @classproperty
        def definite_chars(cls):
            return set([hashes[k] for k in hashes])


        @classproperty
        def default_gap_char(cls):
            return '-'

        @classproperty
        def gap_chars(cls):
             return set('-.')
         

    target_obj = CustomSequence(coded_target)
    cluster_obj = CustomSequence(coded_cluster)
    
    substitution_matrix = make_identity_substitution_matrix(
        match_score=1,
        mismatch_score=-1,
        alphabet=letters
    )
    
    
    alignment = global_pairwise_align(
            target_obj, 
            cluster_obj,
            gap_open_penalty=1,
            gap_extend_penalty=1,
            substitution_matrix=substitution_matrix
        )
    
    return alignment



























    
