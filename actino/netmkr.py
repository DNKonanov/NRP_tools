#!/bin/env python

import sys
import networkx as nx
import matplotlib.pyplot as plt
from math import log

## Define colors
col = {}
col['Streptomyces'] = '#33CCFF'
col['Micromonospora'] = '#6666FF'
col['Salinispora'] = '#00CC00'
col['Mycobacterium'] = '#000099'
col['Gordonia'] = '#9900CC'
col['Rhodococcus'] = '#CC9900'
col['Nocardia'] = '#CC00CC'
col['Amycolatopsis'] = '#CCCC00'
col['base'] = '#333333'

## Read in the edge table
## expects format:  Source Target dist Source_Genus Target_Genus
seen = {}
glist = {}
genusof = {}
for genus in col:
    glist[genus] = []
g = nx.Graph()
f = open(sys.argv[1], "r")
next(f)
for line in f:
    l = line.rstrip().split("\t")
    if(l[0] not in seen.keys()):
        if(l[3] in col.keys()):
            g.add_node(l[0], color=col[l[3]])
            glist[l[3]].append(l[0])
        else:
            g.add_node(l[0], color=col['base'])
            glist['base'].append(l[0])
        seen[l[0]] = 1
        genusof[l[0]] = l[3]
    if(l[1] not in seen.keys()):
        if(l[4] in col.keys()):
            g.add_node(l[1], color=col[l[4]])
            glist[l[4]].append(l[1])
        else:
            g.add_node(l[1], color=col['base'])
            glist['base'].append(l[1])
        seen[l[1]] = 1
        genusof[l[1]] = l[4]
    g.add_edge(l[0], l[1], weight=float(l[2]))
f.close()

layout = 0
if(layout == 1):
    pos=nx.nx_pydot.graphviz_layout(g, prog="neato")
    for genus in glist:
        nx.draw_networkx_nodes(g,pos,node_size=10,nodelist=glist[genus],node_color=col[genus])
    nx.draw_networkx_edges(g,pos,alpha=0.5,width=2)
    nx.write_graphml(g, "actino.graphml")
    plt.axis('off')
    plt.savefig("actino.png", dpi=600)

## Calculate the shannon entropy of each component larger
## than a defined threshold of nodes
thresh = 5
sg = nx.connected_component_subgraphs(g)
c = 1
insubgraphs = 0
for comp in sg:
    if(len(comp) >= thresh):
        gc = {}
        for n in comp.nodes():
            sys.stderr.write("\t".join(["comp"+str(c), n+"\n"]))
            if(genusof[n] in gc.keys()):
                gc[genusof[n]] += 1
            else:
                gc[genusof[n]] = 1
        shannon = 0
        major = {}
        for g in gc:
            if('count' in major.keys()):
                if(major['count'] < gc[g]):
                    major['genus'] = g
                    major['count'] = gc[g]
            else:
                major['genus'] = g
                major['count'] = gc[g]
            prop = float(gc[g]) / len(comp)
            shannon += float(prop*log(prop, 2))
        if(shannon != 0):
            shannon = -shannon
        print "\t".join([str(len(comp)), str(shannon), major['genus']])
        c += 1
    insubgraphs += len(comp)

sys.stderr.write(str(insubgraphs)+" nodes (BGCs) represented in network\n")
 
