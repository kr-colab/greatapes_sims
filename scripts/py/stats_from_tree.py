'''
input: tree_path(str), spp(str), rand_id(str), win_size(int), L(int)
output: csv with all stats(long) for all windows (append existing)
'''
import sys
from stats_funcs import *

def single_pop_stats_from_ts(ts_path, L, win_size):
    #getting the identifier of the treeseq
    #getting all pairwise combinations of pops
    print("entrei")
    ts = pyslim.load(ts_path).simplify()
    s1 = timer()
    print("vou pegar os acs")
    ac, pos = acs_from_ts(ts, 1)
    tmp = pd.DataFrame()
    print("Calculating single population stats...", flush=True)
    pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=win_size,   start=1, stop=L)
    D, windows, counts = allel.windowed_tajima_d(pos, ac, size=win_size, start=1, stop=L)
    tmp['start'] = windows[:,0]
    tmp['end'] = windows[:,1]
    tmp['n_acc'] = n_bases
    tmp['pi_p'+str(j)] = pi
    tmp['tajd_p'+str(j)] = D
    s2 = timer()
    print(("Calculated single pop stats... Time elapsed (min):"+str(round((s2-s1)/60,    3))), flush=True)
    s1=timer()
    return(tmp)

print("input: tree_path(str), spp(str), rand_id(str), win_size(int), L(int)")
assert len(sys.argv) == 6, "Not enough input was provided."

ts_path = sys.argv[1]
spp = sys.argv[2]
rand_id = sys.argv[3]
win_size = sys.argv[4]
L = sys.argv[5]

