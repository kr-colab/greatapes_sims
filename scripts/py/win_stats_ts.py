import allel
import subprocess, msprime, pyslim
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import re
from glob import glob
import pickle
import sys
import itertools
from timeit import default_timer as timer
import pandas as pd

def acs_from_ts(ts, n_pops, N):
    '''
    This function takes a tree sequence,  and returns tuple with a list of allele counts for each subpop and the positions'''
    acs=[]
    hap = allel.HaplotypeArray(ts.genotype_matrix())
    geno = hap.to_genotypes(ploidy=2)
    for i in range(n_pops):
        subpop_indexes = list(np.arange(i*N,(i+1)*N))
        acs.append(geno.count_alleles(subpop=subpop_indexes))
    pos=np.array([s.position for s in ts.sites()])
    return(acs, pos)

def get_meta(ts_path, meta_path):
    matches = re.match( r'.+RAND_(.+).trees', ts_path)
    rand_id = matches.groups()[0]
    files = glob(meta_path+rand_id+"_jid_*.info") #there should be only one match
    meta_fname = files[0]
    meta = pd.read_csv(meta_fname, sep="\t")
    N = int(meta['N1'])
    L = int(meta['L'])
    return(rand_id, N, L)

def win_stats_from_ts(ts_path, rand_id, n_pops, N, L, win_size):
    #getting the identifier of the treeseq
    #getting all pairwise combinations of pops
    x = np.arange(n_pops)
    combs = list(itertools.combinations(x, 2))
    ts = pyslim.load(ts_path).simplify()
    s1 = timer()
    acs, pos = acs_from_ts(ts, n_pops, N)
    tmp = pd.DataFrame()
    print("Calculating single population stats...", flush=True)
    for j in range(n_pops):
        pi, windows, n_bases, counts = allel.windowed_diversity(pos, acs[j], size=win_size,   start=1, stop=L)
        D, windows, counts = allel.windowed_tajima_d(pos, acs[j], size=win_size, start=1,     stop=L)
        if (tmp.empty):
            tmp['start'] = windows[:,0]
            tmp['end'] = windows[:,1]
            tmp['n_acc'] = n_bases
        tmp['pi_p'+str(j)] = pi
        tmp['tajd_p'+str(j)] = D
    s2 = timer()
    print(("Calculated single pop stats... Time elapsed (min):"+str(round((s2-s1)/60,    3))), flush=True)
    s1=timer()
    for k in range(len(combs)):
        dxy, windows, n_bases, counts = allel.windowed_divergence(pos, acs[combs[k][0]],      acs[combs[k][1]], size=win_size, start=1, stop=L)
        tmp['dxy_p'+str(combs[k][0])+'_'+str(combs[k][1])] = dxy
        #fstat, windows, counts = allel.windowed_hudson_fst(pos, acs[combs[k][0]],             acs[combs[k][1]], size=win_size, start=1, stop=L)
    s2 = timer()
    print(("Calculated windowed Dxy... Time elapsed (min):"+str(round((s2-s1)/60,    3))), flush=True)
    return(tmp)

out_path = "/home/murillor/projects/greatapes_sims/output/"
meta_path = "/home/murillor/projects/greatapes_sims/meta/sims/"
prefix = "neut_split_mut"

tree_files = glob(out_path+prefix+"_RAND_*.trees")
for ts_path in tree_files:
    N, L = get_meta(ts_path, meta_path)
    cmd = "sbatch "+rand+".sh"
    out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    if (out.returncode != 0):
        print("sh failed")
        print(out)


sys.stdout.flush()


