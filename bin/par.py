#!/bin/env python

import multiprocessing
import sys, os, re, os.path

def sandpuma_mt(fn,bd):
    newfn = os.path.splitext(fn)[0]
    pref = newfn.split('/')[-1]
    newfn = newfn+os.path.splitext(fn)[1]
    #print 'SANDPUMA parallel instance: '+pref
    ## Make subdir
    if(os.path.isdir(pref)==False):
        os.system('mkdir '+pref)
    ## Copy subset to the subdir
    initcp = ' '.join(['cp', fn, pref])
    os.system(initcp)
    ## Run SANDPUMA
    apcmd = 'perl '+bd+'/allpred_nodep_par.pl '+pref+'/'+newfn
    os.system(apcmd)
    
jobs = []
prefs = []
b = './'
for i,a in enumerate(sys.argv):
    if(i==0):
        b = os.path.dirname(os.path.realpath(a))
        continue
    newfn = os.path.splitext(a)[0]
    pref = newfn.split('/')[-1]
    p = multiprocessing.Process(target=sandpuma_mt, args=(a,b))
    jobs.append(p)
    prefs.append(pref)
    p.start()

#for p in prefs:
#    print p
