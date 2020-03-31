'''
input: tree_path(str), out_path(str), spp(str), rand_id(str), win_size(int), L(int), n(int), center(bool)
path to the tree sequence
path to the output where single_pop_all.tsv will be saved
spp is the species name to annotate the table
rand_id the identifier to the simulation params
rep replicate number
win_size is the window size
L is the chromosome size
n is the size of the sample with which stats will be calculated
center determines if only one centered window will be used to calculate stats on
output: csv with all stats(long) for all windows (append existing)
'''
import sys
import os.path
from stats_funcs import *
import numpy as np
print("input: tree_path(str), filename(str), spp(str), rand_id(str), rep(str), win_size(int), L(int), n(int), center(bool)")
assert len(sys.argv) < 9, "Not enough input was provided."
center=False
ts_path = sys.argv[1]
filename = sys.argv[2]
spp = sys.argv[3]
rand_id = sys.argv[4]
rep = sys.argv[5]
win_size = int(sys.argv[6])
L = int(sys.argv[7])
n = int(sys.argv[8])
if len(sys.argv > 9:
    center = True
print(rand_id)
print(n)
stats = single_pop_stats_from_ts(ts_path, L, win_size, n, center)
stats['spp'] = spp
stats['rand_id'] = rand_id
stats['rep'] = rep
stats.to_csv("single_pop_all.tsv", header=(not os.path.exists("single_pop_all.tsv")), index=False, mode="a")
stats.to_csv(filename, sep="\t", header=(not os.path.exists(filename)), index=False,         mode="a")
