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

def remove_mutations(ts, start, end, proportion):
    '''
    This function will return a new tree sequence the same as the input,
    but after removing each non-SLiM mutation within regions specified in lists
    start and end with probability `proportion`, independently. So then, if we
:/name
    '''
    new_tables = ts.dump_tables()
    new_tables.mutations.clear()
    mutation_map = [-1 for _ in range(ts.num_mutations)]
    for j, mut in enumerate(ts.mutations()):
        keep_mutation = True
        for i in range(len(start)):
            left = start[i]
            right = end[i]
            assert(left < right)
            if i > 0:
                assert(end[i - 1] <= left)
            if mut.position >= left and mut.position < right and len(mut.metadata) == 0:
                keep_mutation = (random.uniform(0, 1) > proportion)
        if keep_mutation:
            mutation_map[j] = new_tables.mutations.num_rows
            if mut.parent < 0:
                new_parent = -1
            else:
                new_parent = mutation_map[mut.parent]
            new_tables.mutations.add_row(site = mut.site, node = mut.node,
                    derived_state = mut.derived_state,
                    parent = new_parent,
                    metadata = mut.metadata)
    return new_tables.tree_sequence()

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
    ts_mut.dump(ts_path_mut)
    s2 = timer()
    print(("Dumped overlaid trees to file", ts_path_mut, "... Time elapsed (min):"+str(round((s2-s1)/60,3))), flush=True)

ts_path = sys.argv[1]
ts_path_mut = sys.argv[2]
neut_mut = sys.argv[3]

overlay_varmut(ts_path, ts_path_mut, neut_mut)
