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

def remove_mutations(ts, starts, ends, prop):
    '''
    This function will return a new tree sequence the same as the input,
    but after removing each non-SLiM mutation within regions specified in 
    lists of start and end positions with probability `proportion`, independently. 
    So then, if we want to add neutral mutations with rate 1.0e-8 within the regions 
    and 0.7e-8 outside the regions, we could do
      ts = pyslim.load("my.trees")
      first_mut_ts = msprime.mutate(ts, rate=1e-8)
      mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)
    :param float proportion: The proportion of mutations to remove.
    '''
    pos = ts.tables.sites.position #getting the positions of all sites
    is_msp = (np.diff(ts.tables.mutations.metadata_offset) == 0) #getting which mutations are from msprime
    #but we want to know which sites are from msprime
    is_msp_site = np.repeat(False, ts.num_sites)
    is_msp_site[ts.tables.mutations.site] = is_msp
    #finding which sites are inside the regions
    breaks=np.concatenate(([-1], starts, ends))
    breaks.sort()
    #np.search sorted is going to return even numbers if the pos is inside one of the regions
    in_regions = np.searchsorted(breaks,pos,"right")%2 == 0
    removable_sites = np.where(np.logical_and(in_regions, is_msp_site))[0]
    #find sites to remove with probability prop
    remove = np.where(np.random.binomial(1,prop,len(removable_sites))==1)[0]
    new_table = ts.tables
    new_table.delete_sites(remove)
    return(new_table.tree_sequence())

def overlay_varmut(ts_path, ts_path_mut, neut_mut, intervals = False):
    s1=timer()
    ts_slim = pyslim.load(ts_path).simplify()
    ts_mut = msprime.mutate(ts_slim, neut_mut, keep=True)
    print("Mutated", ts_path, "in msprime...", flush=True)
    #start, end, del_mut = extract_meta(fpath, gene_info, type)
    '''if del_mut > 0:
        s1=timer()
        ts = remove_mutations(ts_mut, start, end, del_mut/neut_mut)
        s2 = timer()
        print(("Removed extra mutations in genic regions from", fpath, "... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)
    else:
        ts = ts_mut
    s1=timer()'''
    print(ts_path_mut)
    ts_mut.dump(ts_path_mut)
    s2 = timer()
    print(("Dumped overlaid trees to file", ts_path_mut, "... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

ts_path = sys.argv[1]
ts_path_mut = sys.argv[2]
neut_mut = float(sys.argv[3])


overlay_varmut(ts_path, ts_path_mut, neut_mut)
