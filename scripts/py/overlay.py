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



def overlay_varmut(in_ts_path, out_ts_path, mut_rate, recapN, rec_hap_path, ex_file_path, intervals = False):
    s1=timer()
    ts_slim = pyslim.load(in_ts_path)
    if recapN > 0:
        recomb_map = msprime.RecombinationMap.read_hapmap(rec_hap_path)
        ts_slim = ts_slim.recapitate(recombination_map=recomb_map, Ne=recapN)
        print("Recapitated", ts_path, "with pyslim...", flush=True)
    ts_mut = msprime.mutate(ts_slim, mut_rate, keep=True)
    if ex_file_path:
    print("Mutated", in_ts_path, "in msprime...", flush=True)
    ts_mut.dump(out_ts_path)
    s2 = timer()
    recap_message = "" if recap == "" else " and recapped"
    print(("Dumped overlaid"+recap_message+" trees to file", ts_path_mut, "... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

in_ts_path = sys.argv[1]
out_ts_path = sys.argv[2]
mut_rate = float(sys.argv[3])
recapN = int(sys.argv[4])
rec_hap_path = sys.argv[5]
ex_file_path = sys.argv[6]

overlay_varmut(in_ts_path, out_ts_path, mut_rate, recapN, rec_hap_path)
