#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import tskit
import msprime
import numpy as np
import argparse


# In[ ]:


parser = argparse.ArgumentParser(description='Gets stats from unioned tree sequence')
parser.add_argument('infilepath', type=str)
parser.add_argument('outfilepath', type=str)
parser.add_argument('win_size', type=lambda x: int(float(x)))
parser.add_argument('total_mut_rate', type=float)
parser.add_argument('sigma', type=float, help="rates are sampled from (1+N(0,sigma))*total_mut_rate")
parser.add_argument('coords_dict', type=str, help="String of a dictionary with padded and non-padded start and ends of the chromosomic region. Assumes one chromosome only!")
parser.add_argument('--seed', type=int, default=2735, required=False)


# In[ ]:


"""
args = vars(parser.parse_args(["/home/murillor/projects/greatapes_sims/output/4QN9GACHUG2OXET/4QN9GACHUG2OXET_rep0.union.recap.mut.trees",
                               "/home/murillor/projects/greatapes_sims/output/4QN9GACHUG2OXET/4QN9GACHUG2OXET_rep0.union.recap.mut.trees2",
                               "1000000",
                               "2e-8",
                               "0.05",
                               "{'chr':'chr12', 'start':40000000, 'end':50000000, 'padded_start':40000000, 'padded_end': 50000000}"]))
"""


# In[ ]:


args = vars(parser.parse_args())


# In[ ]:


coords= eval(args['coords_dict'])


# In[ ]:


start = coords["start"] - coords["padded_start"]
end = coords["end"]-coords["start"] + start
clen = coords["padded_end"]-coords["padded_start"]
print(end, start)


# In[ ]:


# making sure we can fit the windows into the chunk
# this would fail if the simulated chunk (excluding padding) is not a multiple of the win_size


# In[ ]:


breakpoints = []
if start > 0:
    breakpoints.append(0)
breakpoints += list(range(start,clen+1, args["win_size"]))
assert breakpoints[-1] <= clen
if breakpoints[-1] < clen:
    breakpoints.append(clen)
print(breakpoints)


# In[ ]:


rng = np.random.default_rng(args["seed"])


# In[ ]:


coefs = rng.normal(0,args["sigma"], len(breakpoints)-1)
rates = (1+coefs)*args["total_mut_rate"]


# In[ ]:


rmap = msprime.RateMap(position=breakpoints, rate=rates)


# In[ ]:


ts = tskit.load(args["infilepath"])


# In[ ]:


model = msprime.SLiMMutationModel(type=1, next_id=1)


# In[ ]:


# setting keep=False will erase any existing mutation from the tree sequence.
tsmut = msprime.sim_mutations(ts, model=model, rate=rmap, keep=False)


# In[ ]:


tsmut.dump(args["outfilepath"])

