#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
#%load_ext memory_profiler


# In[2]:


parser = argparse.ArgumentParser(description='Gets stats from unioned tree sequence')
parser.add_argument('infile', type=str)
parser.add_argument('popsfile', type=str)
parser.add_argument('outfile', type=str)
parser.add_argument('sample_size', type=int)
parser.add_argument('coords_dict', type=str, help="String of a dictionary with padded and non-padded start and ends of the chromosomic region. Assumes one chromosome only!")
parser.add_argument('--seed', type=int, default=8991, required=False)


# In[3]:


#%%memit
#args = vars(parser.parse_args(["../../output/GDVOP9EMEV6ZLEG/GDVOP9EMEV6ZLEG_rep0.union.recap.mut.trees","../../output/GDVOP9EMEV6ZLEG/GDVOP9EMEV6ZLEG_rep0.pops","test.txt","10","{'chr': 'chr12', 'start': 60000000, 'end': 70000000, 'padded_start': 60000000, 'padded_end': 70000000}"]))
args = vars(parser.parse_args())
coords_dict = eval(args['coords_dict'])
rng = np.random.default_rng(args['seed'])


# In[5]:


#%%memit
assert os.path.exists(args["infile"]), f"Tree sequence file does not exist {args['infile']}"
assert os.path.exists(args["popsfile"]), f".pops file does not exist {args['popsfile']}"
with open(args["popsfile"], "r") as f:
    pops = eval(f.readline())
recap_tsu = tskit.load(args["infile"])


# In[36]:


#%%memit
# keeping only the focal window
start = coords_dict['start']-coords_dict['padded_start']
stop = (coords_dict['end']-coords_dict['padded_start'])
recap_tsu = recap_tsu.keep_intervals([[start,stop]])


# In[6]:


#%%memit
samples = sample_from_ts(recap_tsu, sample_size=args["sample_size"], rng=rng)


# In[7]:


sample_sets = list(samples.values())


# In[8]:


# getting num samples
n = np.array([len(s) for s in sample_sets], dtype='float')
num_pops = len(sample_sets)

# all two-way comparisons between all pops
twoway = [[x,y] for x in range(num_pops) for y in range(num_pops) if x>=y]
# all four-way comparisons between all pops
fourway = np.array([(twoway[xx] + twoway[yy]) for xx in range(len(twoway)) for yy in range(len(twoway)) if xx>=yy], dtype='int')
# note each row in fourway contains the four pop indices
i, j, k, l = [fourway[:,x] for x in range(4)]


# In[9]:


# covariance between two PPP random variables is just the integral of their products!
# see https://www.sciencedirect.com/science/article/pii/S0040580918301667
# This function takes x, a vector with counts, uses 5 vars from the global env
# n, with the total counts, ijkl with the corresponding indices of pops 1-4.
# Returns the covariance between all four-way combinations of the pops.
def pidxy_cov(x):
    numer = x[i] * (n[j] - x[j]) * x[k] * (n[l] - x[l])
    denom = n[i] * (n[j] - (i == j)) * n[k] * (n[l] - (k == l))
    return numer / denom


# In[10]:


#%%memit
covs = recap_tsu.sample_count_stat(sample_sets=sample_sets, f=pidxy_cov, output_dim=fourway.shape[0])


# In[11]:


cov_labels = np.array(pops)[fourway]


# In[54]:


covdf = pd.DataFrame(cov_labels, columns=["spp2_1", "spp1_1", "spp2_2", "spp1_2"])


# In[55]:


covdf["cov"] = covs


# In[56]:


covdf = covdf.join(pd.DataFrame(coords_dict, index=covdf.index))


# In[61]:


covdf = covdf[list(coords_dict.keys())+["spp1_1", "spp2_1", "spp1_2","spp2_2", "cov"]]


# In[ ]:


covdf.to_csv(args["outfile"], sep="\t", index=False)


# In[ ]:


"""
import collections


rank_counts = collections.Counter(t.rank() for t in tssimp.trees())

most_common = rank_counts.most_common(100)

most_common[0][0]

trees = [tskit.Tree.unrank(tssimp.num_samples, mc[0]) for mc in most_common]

tssimp.num_trees

from IPython.display import SVG, display, HTML
for i in range(10):
    print(most_common[i])
    display(SVG(trees[i].draw_svg(node_labels = node_labels)))

multi_tree_str = ""

for tree in trees:
    multi_tree_str += tree.as_newick(include_branch_lengths=False, node_labels = node_labels)
    multi_tree_str += "\n"

multi_tree_str

import toytree

mtre0 = toytree.mtree(multi_tree_str)

canvas, axes, mark = mtre0.draw_cloud_tree(
    edge_style={
        "stroke-opacity": 0.1,
        "stroke-width": 1,
    },
);
"""

