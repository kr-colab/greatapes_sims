'''
input: tree_path(str), out_path(str), spp(str), rand_id(str), win_size(int), L(int)
output: csv with all stats(long) for all windows (append existing)
'''
import sys
import os.path
from stats_funcs import *
import numpy as np
print("input: tree_path(str), filename(str), spp(str), rand_id(str), rep(str), win_size(int), L(int), n(int), center(bool)")
# if center is True, then stats will be calculated only for one window center in the chr L
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
