#!/bin/env python

from Bio import Phylo
import sys, re

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def get_child_leaves(tree, parent):
    child_leaves = []
    for leaf in tree.get_terminals():
        for node in tree.get_path(leaf):
            if(node == parent):
                child_leaves.append(leaf)
    return child_leaves

def calcscore(scaleto, distance):
    if(distance >= scaleto):
        return 0
    else:
        return float(scaleto - distance) / scaleto

def getscore(scaleto, nd, dist2q, leaf, o):
    score = 0
    nnspec = leaf[o[0]]['spec']
    for n in o:
        curspec = leaf[o[0]]['spec']
        if(nnspec == curspec):
            tmpscore = calcscore(scaleto, float(dist2q[n]) / nd)
            if(tmpscore > 0):
                score += tmpscore
            else:
                break
        else:
            break
    return score

def deeperdive(query, tree, n, p, l):
    ## Want q to nearest dist to be less than nearest to pannearest dist
    q_to_n = tree.distance(l[query]['node'], l[n]['node'])
    n_to_pn = tree.distance(l[n]['node'], l[p]['node'])
    q_to_pn = tree.distance(l[query]['node'], l[p]['node'])
    if((q_to_n < n_to_pn) and (l[n]['spec'] == l[p]['spec'])):
        return (l[n]['spec'], l[query]['id'], 'no_force_needed')
    elif((q_to_n == q_to_pn) and (l[n]['spec'] != l[p]['spec'])):
        return ('no_confident_result', 'NA', 'no_confident_result')
    else:
        parent = get_parent(tree, l[query]['node'])
        sisdist = {}
        for sis in get_child_leaves(tree, parent):
            sisdist[sis.name] = {}
            sisdist[sis.name]['dist'] = tree.distance(l[query]['node'], sis)
            sisdist[sis.name]['node'] = sis
        o = sorted(sisdist, key=sisdist.get)
        for n in o:
            if(n != l[query]['id']):
                fc = re.split("_+", n)
                return ('no_confident_result', 'NA', fc[-1])
        return ('no_confident_result', 'NA', 'no_confident_result')
    
def checkclade(query, lo, hi, wc, tree, l):
    if((lo in l) and (hi in l)): ## Not first or last
        if(l[lo]['spec'] == wc): ## lower bound is wildcard
            checkclade(query, lo-1, hi, wc, l)
        elif(l[hi]['spec'] == wc): ## upper bound is wildcard
            checkclade(query, lo, hi+1, wc, l)
        else:
            ## Get the lca's descendants and check specs
            lca = tree.common_ancestor(l[lo]['node'], l[hi]['node'])
            spec = ''
            iname = ''
            passflag = 1
            for child in get_child_leaves(tree, lca):
                ida = re.split("_+", child.name)
                if(spec != ''):
                    if((ida[-1] != spec) and (ida[-1] != wc)):
                        passflag = 0
                    else:
                        spec = ida[-1]
                        iname = ida[0]
                else:
                    spec = ida[-1]
                    iname = ida[0]
            if(passflag == 0):
                return('deeperdive', 'NA')
            else:
                return(spec, iname)
    else: ## First or last
        return('deeperdive', 'NA')

## set wildcard for query seqs in the tree
wild = 'UNK'
zero_thresh = 0.005

## Set normalization twin pairs
nppref = ['Q70AZ7_A3', 'Q8KLL5_A3']
npnode = ['','']

## Set normalized distance cutoff for nearest neighbor
## (empirically derived default = 2.5)
dcut = 2.5;

## Read tree
treef = sys.argv[1]
tree = Phylo.read(treef, 'newick')

query = []
leaves = {}

## Loop through leaves to find groups
last = ''
group = 1
leafnum = 1

for leaf in tree.get_terminals():
    if(bool(re.search('^'+nppref[0], leaf.name))):
        npnode[0] = leaf
    elif(bool(re.search('^'+nppref[1], leaf.name))):
        npnode[1] = leaf
    ida = re.split("_+", leaf.name)
    if(last != ''): ## Every pass except the first
        if((last != ida[-1]) or (last != wild)): ## begin new group
            group += 1
    last = ida[-1]
    leaves[leafnum] = {}
    leaves[leafnum]['group'] = group
    leaves[leafnum]['id'] = leaf.name
    leaves[leafnum]['spec'] = ida[-1]
    leaves[leafnum]['node'] = leaf
    ## Record queries
    if(ida[-1] == wild):
        query.append(leafnum)
    leafnum += 1
for q in query:
    ## Get distances to knowns
    distfromq = {}
    for n in leaves:
        if((q != n) and (leaves[n]['spec'] != wild)):
            distfromq[n] = tree.distance(leaves[q]['node'], leaves[n]['node'])
    # Sort distances
    o = sorted(distfromq, key=distfromq.get)
    ## Get zero distances
    z = []
    for n in o:
        if(distfromq[n] <= zero_thresh):
            if(distfromq[n] >= 0):
                z.append(n)
        else:
            break
    forcedpred = 'no_force_needed'
    pred = 'no_pred_made'
    hit = 'NA'
    if(len(z) > 0): ## query has zero length known neighbors
        pred = leaves[z[0]]['spec']
        hit = re.search("^(\S+)_.+$", leaves[z[0]]['id']).groups()[0]
    else:
        ## check it out
        pred, hit = checkclade(q, q-1, q+1, wild, tree, leaves)
        if(pred == 'deeperdive'):
            pred, hit, forcedpred = deeperdive(q, tree, o[0], o[1], leaves)
            if(hit != 'NA'):
                hit = re.search("^(\S+)_.+$", hit).groups()[0]
    normdist = tree.distance(npnode[0], npnode[1])
    nn = float(distfromq[o[0]]) / normdist
    nnscore = 0
    snnscore = 0
    if(nn < dcut):
        snnscore = getscore(dcut, normdist, distfromq, leaves, o)
        nnscore = calcscore(dcut, nn)
    print "\t".join([leaves[q]['id'], pred, hit, forcedpred, str(nn), str(nnscore), str(snnscore)])
