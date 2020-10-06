#!/bin/env python

import sys
import networkx as nx
import matplotlib.pyplot as plt
from math import log

seen = {}
g = nx.Graph()
f = open(sys.argv[1], "r")
next(f)
for line in f:
    l = line.rstrip().split("\t")
    if(l[0] not in seen.keys()):
        g.add_node(l[0])
        seen[l[0]] = 1
    if(l[1] not in seen.keys()):
        g.add_node(l[1])
        seen[l[1]] = 1
    g.add_edge(l[0], l[1])
f.close()

sg = nx.connected_component_subgraphs(g)
c = 1
for comp in sg:
    print "\t".join([str(c), str(len(comp))])
    c += 1
