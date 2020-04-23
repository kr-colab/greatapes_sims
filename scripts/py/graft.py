import sys
import os.path
import numpy as np
import pyslim
from stats_funcs import *

ts1=pyslim.load("../../output/eastern-chimp_ELK3ZHJZS3WAM20_rep0.trees")
ts2=pyslim.load("../../output/western-chimp_ELK3ZHJZS3WAM20_rep0.trees")

def get_slim_ids(ts):
    return np.array([n.metadata.slim_id for n in ts.nodes()])

def get_slim_gens(ts):
    return np.array([p.slim_generation for p in ts.slim_provenances])

# getting SLiM generations from the provenances
slim_gens1=get_slim_gens(ts1)
slim_gens2=get_slim_gens(ts2)
assert (len(slim_gens1) == len(slim_gens2)) and (slim_gens1[-1]==slim_gens2[-1]) and (slim_gens1[0]==slim_gens2[0]),
'Based on SLiM generations from the provenances, it does not seem these tree sequences have shared history'

# finding where histories diverge
is_eq = slim_gens1==slim_gens2
first_diff = np.where(np.diff(is_eq) > 0)[0][0]
split_time = abs(slim_gens1[first_diff]-slim_gens2[-1])

# first, identifying equivalency between nodes of the two treeseqs
slim_ids1= get_slim_ids(ts1)
slim_ids2= get_slim_ids(ts2)

n_times2 = np.array([n.time for n in ts2.nodes()])

