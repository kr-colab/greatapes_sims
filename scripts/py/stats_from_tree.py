'''
input: tree_path(str), out_path(str), spp(str), rand_id(str), win_size(int), L(int)
output: csv with all stats(long) for all windows (append existing)
'''
import sys
import os.path
from stats_funcs import *

print("input: tree_path(str), out_path(str), spp(str), rand_id(str), win_size(int), L(int), n(int)")
assert len(sys.argv) == 8, "Not enough input was provided."

ts_path = sys.argv[1]
out_path = sys.argv[2]
spp = sys.argv[3]
rand_id = sys.argv[4]
win_size = int(sys.argv[5])
L = int(sys.argv[6])
n = int(sys.argv[7])
stats = single_pop_stats_from_ts(ts_path, L, win_size, n)
stats['spp'] = spp
stats['rand_id'] = rand_id
filename=out_path+spp+"_single_stats_"+rand_id+".tsv"
stats.to_csv(filename, sep="\t", header=(not os.path.exists(filename)), index=False, mode="a")


