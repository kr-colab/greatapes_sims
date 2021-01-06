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

parser = argparse.ArgumentParser(description='Gets stats from unioned tree sequence')
parser.add_argument('rand_id', type=str)
parser.add_argument('rep', type=str)
parser.add_argument('win_size', type=lambda x: int(float(x)))
parser.add_argument('sample_size', type=int)
parser.add_argument('coords_dict', type=str, help="String of a dictionary with padded and non-padded start and ends of the chromosomic region")
parser.add_argument('--seed', type=int, default=8991, required=False)

args = vars(parser.parse_args())
coords_dict = eval(args['coords_dict'])
out_path = "../../output/"

# Loading tree sequence and list with populations
recap_mut_path = f"{out_path}{args['rand_id']}_rep{args['rep']}.union.recap.mut.trees"
pops_path = f"{out_path}{args['rand_id']}_rep{args['rep']}.pops"
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
contemp_samples = [rng.choice(pop_samples[pid][recap_tsu.tables.nodes.time[pop_samples[pid]] == contemp_time[pid]], args["sample_size"], replace=False) for pid in range(len(pop_samples))]


# windowing
start = coords_dict['start']-coords_dict['padded_start']
stop = (coords_dict['end']-coords_dict['padded_start'])
expected_length = (coords_dict['padded_end']-coords_dict['padded_start'])
assert abs(recap_tsu.sequence_length - expected_length) < 1e-14
# coordinate windows: windows in actual chromosome scale, this only matters if the simulations are not chromosome scale, but instead were windowed
coord_windows = np.arange(start=coords_dict['start'], stop=coords_dict['end']+1, step=args['win_size'])
# windows but in simulation size scale (0-L)
windows = np.arange(start=start,stop=stop+1, step=args["win_size"])
assert coord_windows.shape == windows.shape

# dealing with padding
if start > 0:
    windows = np.insert(windows, 0, [0])
if expected_length > stop:
    windows = np.append(windows, [expected_length])
print(recap_tsu.sequence_length)
print(start, stop)
print(windows, coord_windows)
# obtaining indexes for all possible pairs (including diversity, i.e. i==j for (i,j))
indexes = [(x, y) for x in range(len(pops)) for y in range(len(pops)) if x >= y]

#calculating dxy and diversity
dxy = recap_tsu.divergence(sample_sets=contemp_samples, mode="site", windows=windows, indexes=indexes)
# half matrix + diagonal
assert dxy.shape[1] == ((len(pops)**2 - len(pops))/2) + len(pops)

# getting the species labels
labels = np.array([[pops[i],pops[j]] for i, j in indexes])

# slicing out the padded windows
if start > 0:
    dxy = dxy[1:]
    windows = np.delete(windows, 1)
if expected_length > stop:
    dxy = dxy[:-1]
    windows = windows [:-1]

print(dxy.shape, windows.shape, coord_windows.shape)
assert (dxy.shape[0]+1) == windows.shape[0] == coord_windows.shape[0]

# saving to output
np.savez(f"{out_path}rand-id_{args['rand_id']}_rep_{args['rep']}_win-size_{args['win_size']}_sample-size_{args['sample_size']}.npz", windows=windows, coord_windows=coord_windows, dxy=dxy, labels=labels)
