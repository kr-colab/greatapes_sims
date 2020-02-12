import sys
import os.path
import numpy as np
import pyslim
from stats_funcs import *

print("input: tree_paths_string(str), filename(str), spps_string(str), rand_id(str), rep(str), win_size(int), L(int), n(int)")
assert len(sys.argv) == 9, "More arguments needed than "+str(len(sys.argv))

tree_paths_string = sys.argv[1]
filename = sys.argv[2]
spps_string = sys.argv[3]
rand_id = sys.argv[4]
rep = sys.argv[5]
win_size = int(sys.argv[6])
L = int(sys.argv[7])
n = int(sys.argv[8])

paths = tree_paths_string.split(",")
spps = spps_string.split(",")

all_ac=[]
all_pos=[]
for p in paths:
	ts = pyslim.load(p)
	acs, pos = acs_from_ts(ts, n)
	all_ac.append(acs[0])
	all_pos.append(pos)

x = np.arange(len(all_ac))
combs = list(itertools.combinations(x, 2))
stats = pd.DataFrame()
for k in range(len(combs)):
    tmp = pd.DataFrame()
    ac1=all_ac[combs[k][0]]
    pos1=all_pos[combs[k][0]]
    ac2=all_ac[combs[k][1]]
    pos2=all_pos[combs[k][1]]
    pos = np.unique(np.concatenate((pos1,pos2)))
    pos.sort()
    new_ac1 = np.full((len(pos),2),0)
    new_ac2 = np.full((len(pos),2),0)
    new_ac1[np.where(np.isin(pos,pos1))] = ac1
    new_ac2[np.where(np.isin(pos,pos2))] = ac2
    n_sampled1 = np.sum(ac1,axis=1)[0]
    n_sampled2 = np.sum(ac2,axis=1)[0]
# filling out the empty genotypes, which should be fixed to the reference alleles
# if a different allele fixed, then the genotype matrix return by tskit would be e.g. [1,...,1]
    id_invar1 = np.where(np.sum(new_ac1,axis=1)==0)[0]
    id_invar2 = np.where(np.sum(new_ac1,axis=1)==0)[0]
    new_ac1[id_invar1,0]=n_sampled1
    new_ac2[id_invar2,0]=n_sampled2
    dxy, windows, n_bases, counts = allel.windowed_divergence(pos, new_ac1,  new_ac2, size=win_size, start=1, stop=L)
    tmp['start'] = windows[:,0]
    tmp['end'] = windows[:,1]
    tmp['n_acc'] = n_bases
    tmp['dxy'] = dxy
    tmp['spp1'] = spps[combs[k][0]]
    tmp['spp2'] = spps[combs[k][1]]
    tmp['rand_id'] = rand_id
    tmp['rep'] = rep
    stats=pd.concat([stats,tmp])
stats.to_csv("multi_pop_all.tsv", header=(not os.path.exists("multi_pop_all.tsv")), index=False, mode="a")
stats.to_csv(filename, sep="\t", header=(not os.path.exists(filename)), index=False,         mode="a")
