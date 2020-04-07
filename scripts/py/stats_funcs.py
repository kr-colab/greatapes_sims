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

def acs_from_ts(ts, n=0):
    '''
    This function takes a tree sequence, and returns tuple with a list of allele counts for each subpop, and an array of the positions
    n = sample size
    '''
    #np.random.seed(102849)
    acs=[]
    # sampling n individuals from treeseq
    if(n>0):
        # make sure individuals were alive at the end of sim
        alive = ts.individuals_alive_at(0)
        samp_indivs = np.random.choice(alive, size=n, replace=False)
        samp_nodes = []
        for i in samp_indivs:
            samp_nodes.extend(ts.individual(i).nodes)
        samp_nodes = np.array(samp_nodes)
        print(samp_nodes.shape, "shape samp nodes")
        ts = ts.simplify(samp_nodes)
    print("entrei nos acs")
    print(ts.genotype_matrix())
    hap = allel.HaplotypeArray(ts.genotype_matrix())
    geno = hap.to_genotypes(ploidy=2)
    print("fiz hap and geno matrix")
    # getting the pop for each individual and counting alleles for each pop
    ind_pops = np.array([ts.individual(i).population for i in ts.individuals_alive_at(0)])
    print("sample size", n)
    print(ind_pops.shape, "shape")
    for i in range(len(np.unique(ind_pops))):
        subpop_indexes = np.where(ind_pops==i)[0]
        print(subpop_indexes)
        print(geno.shape)
        acs.append(geno.count_alleles(subpop=subpop_indexes))
    print("ja separei os acs per subpop")
    pos=np.array([s.position for s in ts.sites()])
    print(len(acs))
    return(acs, pos)

def single_pop_stats_from_ts(ts_path, L, win_size, n, pad=0):
    print("entrei")
    ts = pyslim.load(ts_path)
    s1 = timer()
    print("vou pegar os acs")
    acs, pos = acs_from_ts(ts, n)
    ac = acs[0]
    tmp = pd.DataFrame()
    print("Calculating single population stats...", flush=True)

    if pad>0:
        # dealing with the 1st window in chr case
        if L < win_size + 2*pad: {
            start = 1
        } else {
            start = pad + 1
        }
        stop = start + win_size
    else:
        start = 1
        stop = L

    pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=win_size, start=start, stop=stop)
    D, windows, counts = allel.windowed_tajima_d(pos, ac, size=win_size, start=start, stop=stop)
    tmp['start'] = windows[:,0]
    tmp['end'] = windows[:,1]
    tmp['n_acc'] = n_bases
    tmp['pi'] = pi
    tmp['tajd'] = D
    s2 = timer()
    print(("Calculated single pop stats... Time elapsed (min):"+str(round((s2-s1)/60, 3))), flush=True)
    s1=timer()
    return(tmp)

def get_meta(ts_path, meta_path):
    '''
    This function gets the genome sequence length and the random identifier based on a path to a tree file (ts_path) and the path to the metadata folder.
    '''
    matches = re.match( r'.+RAND_(.+).trees', ts_path)
    rand_id = matches.groups()[0]
    meta_fname = meta_path+rand_id+".meta"
    meta = pd.read_csv(meta_fname, sep="\t")
    L = int(meta['L'])
    return(rand_id, L)

def win_stats_from_ts(ts_path, rand_id, n_pops, L, win_size):
    #getting the identifier of the treeseq
    #getting all pairwise combinations of pops
    print("entrei")
    x = np.arange(n_pops)
    combs = list(itertools.combinations(x, 2))
    ts = pyslim.load(ts_path).simplify()
    s1 = timer()
    print("vou pegar os acs")
    acs, pos = acs_from_ts(ts, n_pops)
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

def write_sh(out_path, meta_path, script_path, ts_path, win_size, n_pops, prefix, time = "48:00:00", mem = "64G"):
    matches = re.match( r'.+RAND_(.+).trees', ts_path)
    rand_id = matches.groups()[0]
    sh_name = rand_id+"_"+str(win_size)+".sh"
    with open(sh_name, "w") as fh:
        print("#!/bin/bash", file=fh)
        #SBATCH env variables
        print("#SBATCH --account=kernlab\n#SBATCH --partition=kern\n#SBATCH --job-name="+prefix+"\n#SBATCH --time="+time+"\n#SBATCH --mem "+mem+"\n#SBATCH --open-mode=append"+"\n#SBATCH --output="+meta_path+rand_id+".info"+"\n#SBATCH --error="+meta_path+rand_id+".info", file=fh)
        #modules to load on talapas
        print("\nmodule use /projects/apps/shared/modulefiles/\nmodule load python3 tskit SLiM\n", file=fh)
        #slim command
        print("python "+script_path+" "+out_path+" "+meta_path+" "+ts_path+" "+str(win_size)+" "+str(n_pops)+" "+prefix, file=fh)
        print("seff $SLURM_JOBID", file=fh)
