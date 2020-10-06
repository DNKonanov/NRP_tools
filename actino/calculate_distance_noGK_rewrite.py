#!/usr/bin/env python
#-- Bioinformatics Group @ Wageningen University  --
#-- Adapted 2016 Marc Chevrette
#-- Genetics @ UW-Madison
#-- This is the module to calculate the distance between each CLP pair

import sys, os
from munkres import Munkres
import math
import itertools
import pandas as pd
import numpy as np

## Set params
global Jaccardw
global GKw
global DDSw
Jaccardw = 0.667
DDSw = 0.333
seq_simScore = "./full_seqsim_matrix.txt"
cluster = "./corenrps.faa"
pseudo_aa = "./corenrps.clade.txt"
outfile = "mc_dist.txt"
tree_outfile = "mc_tree_upgmma_all_clade.nwk"

def cluster_distance(A, B, nbhood, f):
    #key is name of GC, values is list of specific pfam domain names
    clusterA = BGCs[A] #will contain a dictionary where keys are pfam domains, and values are domains that map to a specific sequence in the DMS variable
    clusterB = BGCs[B] #{ 'general_domain_name_x' : ['specific_domain_name_1', 'specific_domain_name_2'], etc }
    #A and B are lists of pfam domains

    try:
        #calculates the intersect
        Jaccard = len(set(clusterA.keys()) & set(clusterB.keys())) / \
              float( len(set(clusterA.keys())) + len(set(clusterB.keys())) \
              - len(set(clusterA.keys()) & set(clusterB.keys())))
    except ZeroDivisionError:
        print "Zerodivisionerror during the Jaccard distance calculation. Can only happen when one or more clusters contains no pfam domains."
        print "keys of clusterA", A, clusterA.keys()
        print "keys of clusterB", B, clusterB.keys()

    intersect = set(clusterA.keys() ).intersection(clusterB.keys()) #shared pfam domains
    not_intersect = []
    #DDS: The difference in abundance of the domains per cluster
    #S: Max occurence of each domain

    DDS,S = 0,0
    SumDistance = 0
    pair = ""

    for domain in set(clusterA.keys() + clusterB.keys()):
        if domain not in intersect:
            not_intersect.append(domain)

    for unshared_domain in not_intersect: #no need to look at seq identity or anchors, since these domains are unshared
        #for each occurence of an unshared domain do DDS += count of domain and S += count of domain
        dom_set = []
        try:
            dom_set = clusterA[unshared_domain]
        except KeyError:
            dom_set = clusterB[unshared_domain]

        DDS += len(dom_set)
        S += len(dom_set)


    for shared_domain in intersect:
        seta = clusterA[shared_domain]
        setb = clusterB[shared_domain]
        # print seta, setb

        if len(seta+setb) == 2: #The domain occurs only once in both clusters
            pair = tuple(sorted([seta[0],setb[0]]))

            try:
                SumDistance = 1-DMS[shared_domain][pair][0]
                # print 'SumDistance1', SumDistance

            except KeyError:
                print "KeyError on", pair

                errorhandle = open(str(A)+".txt", 'w')
                errorhandle.write(str(pair)+"\n")
                errorhandle.write(str(DMS[shared_domain]))
                errorhandle.close()
            S += max(len(seta),len(setb))
            DDS += SumDistance

        else:                   #The domain occurs more than once in both clusters
            # accumulated_distance = 0

            DistanceMatrix = [[1 for col in range(len(setb))] for row in range(len(seta))]
            # print DistanceMatrix
            for domsa in range(len(seta)):
                # print seta[domsa]
                for domsb in range(domsa, len(setb)):
                    # print setb[domsb]
                    pair = tuple(sorted([seta[domsa], setb[domsb]]))
                    try:
                        Similarity = DMS[shared_domain][pair][0]
                        # print Similarity
                    except KeyError:
                        print "KeyError on", pair
                        errorhandle = open(str(B)+".txt", 'w')
                        errorhandle.write(str(pair)+"\n")
                        errorhandle.write(str(DMS[shared_domain]))
                        errorhandle.close()
                        Similarity=0

                    seq_dist = 1-Similarity
                    DistanceMatrix[domsa][domsb] = seq_dist


            #Only use the best scoring pairs
            Hungarian = Munkres()
            #print "DistanceMatrix", DistanceMatrix
            BestIndexes = Hungarian.compute(DistanceMatrix)
            # print "BestIndexes", BestIndexes
            accumulated_distance = sum([DistanceMatrix[bi[0]][bi[1]] for bi in BestIndexes])
            # print "accumulated_distance", accumulated_distance
            SumDistance = (abs(len(seta)-len(setb)) + accumulated_distance)  #diff in abundance + sequence distance
            # print 'SumDistance2', SumDistance


            S += max(len(seta),len(setb))
            DDS += SumDistance



    #  calculate the Goodman-Kruskal gamma index
    A_pseudo_seq = list(cluster_seq[A])
    B_pseudo_seq = list(cluster_seq[B])
    Ar = [item for item in A_pseudo_seq]
    Ar.reverse()
    #GK = max([calculate_GK(A_pseudo_seq, B_pseudo_seq, nbhood = nbhood), calculate_GK(Ar, B_pseudo_seq, nbhood = nbhood)])

    DDS /= float(S)
    DDS /= float(S)
    DDS = math.exp(-DDS) #transform from distance to similarity score
    # print 'DDS', DDS



    # DDS = 1-DDS #transform into similarity

    #Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK)
    Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS)
    #Similarity_score = (Jaccardw * Jaccard) + (DDSw * DDS) + (GKw * GK)
    Similarity_score = (Jaccardw * Jaccard) + (DDSw * DDS)
    # Similarity_score = 1 - DDS
    if Distance < 0:
        print "negative distance", Distance, "DDS", DDS, pair
        print "Probably a rounding issue"
        print "Distance is set to 0 for these clusters"
        Distance = 0
    #print A, B, Jaccard, GK, DDS
    f.write("\t".join([str(A), str(B), str(Jaccard), str(DDS)])+"\n")
    return Similarity_score



## Generate BGCs and DMS
input = open(cluster, 'r').read().split('>')[1:]
cluster_id = [i.split('\n')[0] for i in input]
pseudo_seq = [i.split('\n')[1] for i in input]
cluster_seq = dict(zip(cluster_id, pseudo_seq))
print "Generating BGCs"
BGCs = {}
for id in cluster_seq.keys():
    cluster_i = {}
    seq = cluster_seq[id]
    seq_dict = {}
    for i in range(len(seq)):
        seq_dict[id+'_'+str(i+1)] = seq[i]
    domain_name = list(set(seq_dict.values()))
    for d in domain_name:
        key = d
        value = []
        for index, domain in seq_dict.items():
            if domain == d:
                value.append(index)
                cluster_i[key] = value
            BGCs[id] = cluster_i
print "DONE"
DMS = {}
print "Parse seqSim"
DM_input = open(seq_simScore, 'r').read().split('\n')
domain_name = []
for l in DM_input:
    for i in [0,5]:
        try:
            domain_name.append(l.split(",")[0].split("|")[i])
        except:
            pass
domain_index = dict(zip(domain_name, range(0,len(domain_name))))
print "DONE"
print "Generating DMS"
a=0
z=0
aa_list = [a for a in open(pseudo_aa, 'r').read().split('\n')[1:] if a != '']
AA = [m.split('\t')[1] for m in aa_list]
for aa in AA:
    d = os.system("date")
    print str(d)+" "+str(aa)
    domain_w_aa = []
    for i in BGCs.keys():
        cluster_i = BGCs[i]
        cluster_i_aa = cluster_i.keys()
        if aa in cluster_i_aa:
            domain_w_aa = domain_w_aa + cluster_i[aa]
    aa_pair_list = [tuple(sorted(list(i))) for i in list(itertools.combinations(domain_w_aa, 2))]
    aa_pair_dict = {}
    print str(len(aa_pair_list))+" aa pairs"
    if len(aa_pair_list) < 1:
        continue
    for p in aa_pair_list:
        compd1a = p[0].split('_')
        index1 = ['A'+ compd1a.pop()]
        compd2a = p[1].split('_')
        index2 = ['A'+ compd2a.pop()]
        A = ['%s|%s' % t for t in zip([p[0]]*len(index1), index1)]
        B = ['%s|%s' % t for t in zip([p[1]]*len(index2), index2)]
        score = []
        for i in range(len(A)):
            A_i = A[i]
            B_i = B[i]
            row_index = domain_index.get(A_i)
            col_index = domain_index.get(B_i)
            if row_index <= col_index:
                try:
                    row = DM_input[col_index]
                    col = row_index
                    score_list = row.split(",")
                    score.append(1-float(score_list[col+1]))
                except:
                    score = [0]
            else:
                try:
                    row = DM_input[row_index]
                    col = col_index
                    score_list = row.split(",")
                    score.append(1-float(score_list[col+1]))
                except:
                    score = [0]
        score=score[0]
        if score == 0:
            z += 1
        else:
            a += 1
        aa_pair_dict[p] = (score, 0)
    DMS[aa] = aa_pair_dict
print "DONE"
#print "Zero scores: "+str(z)+"\tNonzero scores: "+str(a)
#exit()
print "Calculating distance matrix"
cluster_list = cluster_seq.keys()
scale = 1
nbhood = 3
r = len(cluster_list)
c = len(cluster_list)
Dist = [[0 for x in range(r)] for y in range(c)]
dist_score_assembly_line = pd.DataFrame(Dist, index = cluster_list, columns = cluster_list)
pairs = list(itertools.combinations(cluster_list, 2))
print "there are "+str(len(pairs))+" pairs..."
print "Cluster distance"
f = open("cd_output.txt","w") 
for p in pairs:
    sim_score = cluster_distance(A = p[1],  B=p[0], nbhood=nbhood, f=f)
    rowname = p[1]
    colname = p[0]
    dist_score_assembly_line.ix[rowname, colname] = 1- sim_score/scale
f.close()
print "DONE"
dist_score_assembly_line.to_csv(outfile)

#-- Plot the tree
#import Bio.Phylo
#from Bio.Phylo.TreeConstruction import _Matrix, _DistanceMatrix
#from Bio.Phylo.TreeConstruction import _Matrix, _DistanceMatrix
#from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
#names = cluster_list
#score = [s for s in open(outfile, 'r').read().split('\n')[1:] if s != '']
#matrix = []
#for i in range(len(score)):
#    input_i = score[i].split(',')[1:(i+2)]
#    input_i_int = [float(n) for n in input_i]
#    matrix.append(input_i_int)
#m = _DistanceMatrix(names, matrix)

#constructor = DistanceTreeConstructor()
#tree1 = constructor.upgma(m)
#Bio.Phylo.draw_ascii(tree1)
#Bio.Phylo.write(tree1, tree_outfile, 'newick')


