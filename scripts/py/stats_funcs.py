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

def acs_from_ts(ts, n_pops):
    '''
    This function takes a tree sequence, and returns tuple with a list of allele counts for each subpop, and an array of the positions'''
    acs=[]
    print("entrei nos acs")
    hap = allel.HaplotypeArray(ts.genotype_matrix())
    geno = hap.to_genotypes(ploidy=1)
    print("fiz hap and geno matrix")
    for i in range(n_pops):
        subpop_indexes = ts.samples(population=i).tolist()
        acs.append(geno.count_alleles(subpop=subpop_indexes))
    print("ja separei os acs per subpop")
    pos=np.array([s.position for s in ts.sites()])
    return(acs, pos)

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
