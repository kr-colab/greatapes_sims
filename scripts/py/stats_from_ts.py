import tskit
import pyslim
import msprime
import dendropy
import glob
import numpy as np
import pandas as pd
import warnings
import functools
import argparse
import operator
import os

parser.add_argument('rand_id', type=str)
parser.add_argument('rep', type=str)
parser.add_argument('win_size', type=lambda x: int(float(x)))
parser.add_argument('sample_size', type=int)
parser.add_argument('--seed', type=int, default=8991, required=False)

args = vars(parser.parse_args())

# Loading tree sequence and list with populations
recap_mut_path = f"{trees_path}{args['rand_id']}_rep{args['rep']}.union.recap.mut.trees"
pops_path = f"{trees_path}{args['rand_id']}_rep{args['rep']}.pops"
assert os.path.exists(recap_mut_path) and os.path.exists(pops_path), f"Trees file or .pops file does not exist for {args['rand_id']}_{args['rep']}"
recap_tsu = pyslim.load(recap_mut_path)
with open(pops_path, "r") as f:
    pops = eval(f.readline())

rng = np.random.default_rng(args['seed'])
# getting contemporary samples
# note the time of "contemporary" samples varies bc of differences in generation times
# TODO: sample individuals not nodes
pop_samples = [recap_tsu.samples(population_id=i+1) for i in range(len(pops))]
contemp_time = [np.min(recap_tsu.tables.nodes.time[samples]) for samples in pop_samples]
contemp_samples = [rng.choice(pop_samples[pid][recap_tsu.tables.nodes.time[pop_samples[pid]] == contemp_time[pid]], args["sample_size"], replace=False)
                                                        for pid in range(len(pop_samples))]


# windowing
windows = np.arange(start=0,stop=recap_tsu.sequence_length, step=args["win_size"])
if not np.isclose(recap_tsu.sequence_length, windows[-1], rtol=1e-12):
    windows = np.append(windows, [recap_tsu.sequence_length])

# obtaining indexes for all possible pairs (including diversity, i.e. i==j for (i,j))
indexes = [(x, y) for x in range(len(pops)) for y in range(len(pops)) if x >= y]

#calculating dxy and diversity
dxy = recap_tsu.divergence(sample_sets=contemp_samples, mode="site", windows=windows, indexes=indexes)
# half matrix + diagonal
assert dxy.shape[1] == ((len(pops)**2 - len(pops))/2) + len(pops)

# getting the species labels
labels = np.array([[pops[i],pops[j]] for i, j in indexes])

# saving to output
np.savez(f"{trees_path}rand-id_{args['rand_id']}_rep_{args['rep']}_win-size_{args['win_size']}_sample-size_{args['sample_size']}.npz", windows, dxy, labels)
