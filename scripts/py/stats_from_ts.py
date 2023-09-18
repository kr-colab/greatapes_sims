import tskit
import dendropy
import glob
import numpy as np
import pandas as pd
import warnings
import functools
import argparse
import operator
import os
from helper_functions import *

parser = argparse.ArgumentParser(description='Gets stats from unioned tree sequence')
parser.add_argument('infilepath', type=str)
parser.add_argument('popsfilepath', type=str)
parser.add_argument('outfilepath', type=str)
parser.add_argument('win_size', type=lambda x: int(float(x)))
parser.add_argument('sample_size', type=int)
parser.add_argument('coords_dict', type=str, help="String of a dictionary with padded and non-padded start and ends of the chromosomic region. Assumes one chromosome only!")
parser.add_argument('rand_id', type=str)
parser.add_argument('--delete-exons',help="If passed, all the trees from the exons will be removed", action="store_true")
parser.add_argument('--seed', type=int, default=8991, required=False)

args = vars(parser.parse_args())
coords_dict = eval(args['coords_dict'])

# Loading tree sequence and list with populations
assert os.path.exists(args["infilepath"]) and os.path.exists(args["popsfilepath"]), f"Trees file or .pops file does not exist: {args['infilepath']}, {args['popsfilepath']}"
print("about to load", flush=True)
recap_tsu = tskit.load(args["infilepath"])
with open(args["popsfilepath"], "r") as f:
    pops = eval(f.readline())
print("finished loading", flush=True)
rng = np.random.default_rng(args['seed'])
# getting contemporary samples
# note the time of "contemporary" samples varies bc of differences in generation times
# TODO: sample individuals not nodes
samples = sample_from_ts(recap_tsu, sample_size=args["sample_size"], rng=rng)
sample_sets = list(samples.values())
print("samples selected", flush=True)

if args["delete_exons"]:
    print("Deleting trees from exon intervals", flush=True)
    ex_path = f"../../output/maps/{args['rand_id']}_exons.tsv"
    # exons file
    exons = pd.read_csv(ex_path,sep="\t",header=None)
    # removing extraneous columns?
    exons = exons.iloc[:,:3]
    intervals = exons.iloc[:,1:].to_numpy()
    pos = recap_tsu.tables.sites.position
    is_win_exon = ~(np.searchsorted(intervals.flatten(),pos, side="right") % 2 == 0)
    sites_to_delete = np.where(is_win_exon)[0]
    recap_tsu = recap_tsu.delete_sites(sites_to_delete)

# windowing
start = coords_dict['start']-coords_dict['padded_start']
stop = (coords_dict['end']-coords_dict['padded_start'])
expected_length = (coords_dict['padded_end']-coords_dict['padded_start'])
assert abs(recap_tsu.sequence_length - expected_length) < 1e-14
# coordinate windows: windows in actual chromosome scale, this only matters if the simulations are not chromosome scale, but instead were windowed
coord_windows = np.arange(start=coords_dict['start'], stop=coords_dict['end'], step=args['win_size'])
# windows but in simulation size scale (0-L)
windows = np.arange(start=start,stop=stop, step=args["win_size"])
if coord_windows[-1] != coords_dict['end']:
    coord_windows = np.concatenate([coord_windows, [coords_dict["end"]]])
if windows[-1] != stop:
    windows = np.concatenate([windows, [stop]])
assert coord_windows.shape == windows.shape

# dealing with padding
if start > 0:
    windows = np.insert(windows, 0, [0])
if expected_length > stop:
    windows = np.append(windows, [expected_length])
print(recap_tsu.sequence_length, flush=True)
print(start, stop, flush=True)
print(windows, coord_windows, flush=True)
# obtaining indexes for all possible pairs (including diversity, i.e. i==j for (i,j))
indexes = [(x, y) for x in range(len(pops)) for y in range(len(pops)) if x >= y]

#calculating dxy and diversity
dxy = recap_tsu.divergence(sample_sets=sample_sets, mode="site", windows=windows, indexes=indexes)
# half matrix + diagonal
assert dxy.shape[1] == ((len(pops)**2 - len(pops))/2) + len(pops)

# getting the species labels
labels = np.array([[pops[i],pops[j]] for i, j in indexes])
print(indexes, flush=True)
print(labels, flush=True)
# slicing out the padded windows
if start > 0:
    dxy = dxy[1:]
    windows = np.delete(windows, 1)
if expected_length > stop:
    dxy = dxy[:-1]
    windows = windows [:-1]

# creating array with chrom name
chrom = np.full(shape=dxy.shape[0], fill_value=coords_dict['chr'])
print(dxy.shape, windows.shape, coord_windows.shape, flush=True)
assert (dxy.shape[0]+1) == windows.shape[0] == coord_windows.shape[0]

# saving to output
np.savez(args["outfilepath"], windows=windows[:-1], coord_windows=coord_windows[:-1], chrom=chrom, dxy=dxy, labels=labels)
