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
    Do not simplify your tree before running this function, otherwise you will be
    throwing out more mutations than intended. This is necessary, because you do
    not want to throw out any mutations that happened in the recapitation.
    So then, if we want to add neutral mutations with rate 1.0e-8 within the regions
    and 0.7e-8 outside the regions, we could do
      ts = pyslim.load("my.trees")
      first_mut_ts = msprime.mutate(ts, rate=1e-8)
      mut_ts = remove_mutations(first_mut_ts, start, end, 0.3)
    :param float proportion: The proportion of mutations to remove.
    '''
    pos = ts.tables.sites.position #getting the positions of all sites
    # you don't want to remove mutations that happened in the neutral recapitation period.
    is_post_recap = ts.tables.nodes.time[ts.tables.mutations.node] < ts.slim_generation
    is_msp = (np.diff(ts.tables.mutations.metadata_offset) == 0) #getting which mutations are from msprime
    #but we want to know which sites are from msprime
    is_msp_site = np.repeat(False, ts.num_sites)
    is_msp_site[ts.tables.mutations.site] = is_msp
    #finding which sites are inside the regions
    breaks=np.concatenate(([-1], starts, ends))
    breaks.sort()
    #np.search sorted is going to return even numbers if the pos is inside one of the regions
    in_regions = np.searchsorted(breaks,pos,"right")%2 == 0
    removable_sites = np.where(np.logical_and(np.logical_and(in_regions, is_msp_site), is_post_recap))[0]
    #find sites to remove with probability prop
    remove = np.where(np.random.binomial(1,prop,len(removable_sites))==1)[0]
    new_table = ts.tables
    new_table.delete_sites(remove)
    return(new_table.tree_sequence())

def recap(ts, recapN, rec_hap_path=""):
    if (rec_hap_path==""):
        recapped = ts.recapitate(recombination_rate=1e-8, Ne=int(recapN))
    else:
        recomb_map = msprime.RecombinationMap.read_hapmap(rec_hap_path)
        recapped = ts.recapitate(recombination_map=recomb_map, Ne=int(recapN))
    print("Recapitated with recapN of "+recapN+" using pyslim...", flush=True)
    return(recapped)

def remove_extra(ts, ex_file_path, mut_rate, sel_mut_rate, slim_gen):
    assert slim_gen == ts.slim_generation, "You shouldn't simplify the tree before    running remove_mutations."
    prop_remove = sel_mut_rate/mut_rate
    exons = pd.read_csv(ex_file_path,sep="\t")
    ts_removed = remove_mutations(ts, exons.iloc[:,1].to_numpy(), exons.iloc[:,2])
    print("Removed extra "+str(prop_remove*100)+"% of mutations within the specified      regions...")
    return(ts_removed)

def overlay_varmut(in_ts_path, out_ts_path, mut_rate, recapN="", rec_hap_path="", ex_file_path="", sel_mut_rate=0):
    s1=timer()
    ts_slim = pyslim.load(in_ts_path)
    slim_gen = ts_slim.slim_generation
    if recapN != '':
        ts_slim = recap(ts_slim, recapN, rec_hap_path)
    ts_mut = msprime.mutate(ts_slim, mut_rate, keep=True)
    print("Mutated", in_ts_path, "in msprime...", flush=True)
    if ex_file_path != "" and sel_mut_rate!= 0:
        s2 = timer()
        ts_mut = remove_extra(ts_mut, ex_file_path, mut_rate, sel_mut_rate, slim_gen)
        s3 = timer()
        print("Time elapsed to remove (min):"+str(round((s3-s2)/60,3)), flush=True)
    ts_final = ts_mut.simplify()
    recap_message = "" if recapN == "" else " and recapped"
    s4 = timer()
    ts_final.dump(out_ts_path)
    print(("Dumped overlaid"+recap_message+" trees to file", out_ts_path, "... Time elapsed (min):"+str(round((s4-s1)/60,3))), flush=True)

print("arg1: in_ts_path, arg2: out_ts_path, arg3: mut_rate, arg4:recapN, arg5: rec_hap_path, arg6: ex_file_path, arg7: sel_mut_rate")

assert len(sys.argv) > 3, "Not enough input was provided."

in_ts_path = sys.argv[1]
out_ts_path = sys.argv[2]
mut_rate = float(sys.argv[3])
if (len(sys.argv) > 4):
    recapN = sys.argv[4]
    if (len(sys.argv)>6):
        ex_file_path = sys.argv[6]
        rec_hap_path = sys.argv[5]
        sel_mut_rate = sys.argv[7]

overlay_varmut(in_ts_path, out_ts_path, mut_rate, recapN)#, rec_hap_path), ex_file_path, sel_mut_rate)
