#!/usr/bin/env python
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

# variables
parser = argparse.ArgumentParser(description='Perform the union and output stats from tree sequences.')
parser.add_argument('rand_id', type=str)
parser.add_argument('rep', type=str)
parser.add_argument('total_mut_rate', type=float)
parser.add_argument('region_mut_rate', type=float)
parser.add_argument('sample_size', type=int)
parser.add_argument('win_size', type=lambda x: int(float(x)))
parser.add_argument('--recapN', type=int, default=10000, required=False)
parser.add_argument('--seed', type=int, default=8991, required=False)

args = vars(parser.parse_args())

"""
# rand_id and rep
args = {
    "rand_id": "TZPNGS0UY29NGB3",
    "rep": "0",
    "total_mut_rate": 1e-8,
    "region_mut_rate": 0,
    "recapN": 10000,
    "sample_size": 10,
    "win_size": 10**6,
    "seed": 8297,
}
"""

## metadata paths
rec_hap_path = f"../../meta/maps/{args['rand_id']}_recrate.hapmap"
ex_path = f"../../meta/maps/{args['rand_id']}_exons.tsv"


edges_path = "../../meta/edges_meta.tsv"
sims_sum_path = "../../output/rand_id_params.tsv"
sims_full_path = "../../output/sims_info.tsv"
sims_header_path = "../../output/header_sims_info.tsv"
trees_path = "../../output/"

## loading metadata
# edges contains all the edges and info about N and number of generations
edges = pd.read_csv(edges_path,sep="\t")
edges.parent= edges.parent.fillna("")
edges["edge"] = edges["edge"].str.replace('_','-')
edges["parent"] = edges["parent"].str.replace('_','-')
# sims_sum and sims_full relate args["rand_id"]s to simulation parameters
sims_sum = pd.read_csv(sims_sum_path,sep="\t")
sims_full= pd.read_csv(sims_full_path,sep="\t", header=None)
header = pd.read_csv(sims_header_path,sep="\t")
sims_full.columns = header.columns


# getting all output files and grouping by args["rand_id"] and args["rep"]
tree_files = glob.glob(trees_path+"*[0-9].trees")
pattern = f"{args['rand_id']}_rep{args['rep']}"
n_matches = sum(1 for file in tree_files if pattern in file)
# making sure we got all the files
assert n_matches == edges.shape[0]


# getting the phylo tree adn annotating with branch lengths,
tree = build_tree_from_df(edges)
tree = add_blen_from_meta(tree, sims_full, args["rand_id"])

# performing the union
tsu,  pops = union_tseqs(tree,args["rand_id"],args["rep"])
tsu = pyslim.load_tables(tsu.tables)
print(tsu.slim_generation)
assert tsu.max_root_time.is_integer()
tsu = pyslim.annotate_defaults(tsu, tsu.metadata["SLiM"]["model_type"], int(tsu.max_root_time))
print(tsu.slim_generation)
slim_gen = tsu.slim_generation
# asserting within population coalescen
assert len(set([tsu.node(u).population for t in tsu.trees() for u in t.roots])) == 1
tsu.dump(f"{trees_path}{args['rand_id']}_rep{args['rep']}.union.trees")


recomb_map = msprime.RecombinationMap.read_hapmap(rec_hap_path)
recap_tsu = tsu.recapitate(recombination_map=recomb_map, Ne=args["recapN"])
del tsu # too much ram
print(slim_gen, recap_tsu.max_root_time, recap_tsu.num_mutations)


# In[149]:


mut_map = msp_mutation_rate_map(exons, args["total_mut_rate"], region_mut_rate, int(recap_tsu.sequence_length))
model_recap = msprime.SLiMMutationModel(type=3) # TODO: figure out the type number from the treeseq
model_slim = msprime.SLiMMutationModel(type=4) # TODO: figure out the type number from the treeseq
print("Before mutate:", recap_tsu.num_mutations)
recap_tsu = msprime.mutate(recap_tsu, end_time=slim_gen, model=model_recap, rate=total_mut_rate, keep=True)
print("Mutations added in the recapitation:", recap_tsu.num_mutations)
recap_tsu = msprime.mutate(recap_tsu, start_time=slim_gen, model=model_slim, rate=mut_map, keep=True)
print("Total mutations:", recap_tsu.num_mutations)


rng = np.random.default_rng(args['seed'])
# getting contemporary samples
# note the time of "contemporary" samples varies bc of differences in generation times
pop_samples = [recap_tsu.samples(population_id=i+1) for i in range(len(pops))]
contemp_time = [np.min(recap_tsu.tables.nodes.time[samples]) for samples in pop_samples]
contemp_samples = [rng.choice(pop_samples[pid][recap_tsu.tables.nodes.time[pop_samples[pid]] == contemp_time[pid]], args["sample_size"], replace=False)
                                                        for pid in range(len(pop_samples))]

# windowing
windows = np.arange(start=0,stop=recap_tsu.sequence_length, step=args["win_size"])
if not np.isclose(recap_tsu.sequence_length, windows[-1], rtol=1e-12):
    windows = np.append(windows, [recap_tsu.sequence_length])

# obtaining indexes for all possible pairs (including diversity, i.e. i==j for (i,j))
indexes = [(x, y) for x in range(len(pops)) for y in range(len(pops)) if x >= y]

#calculating dxy and diversity
dxy = recap_tsu.divergence(sample_sets=contemp_samples, mode="site", windows=windows, indexes=indexes)
# half matrix + diagonal
assert dxy.shape[1] == ((len(pops)**2 - len(pops))/2) + len(pops)

# getting the species labels
labels = np.array([[pops[i],pops[j]] for i, j in indexes])

# saving to output
np.savez(f"{trees_path}rand-id_{args['rand_id']}_rep_{args['rep']}_win-size_{args['win_size']}_sample-size_{args['sample_size']}.npz", windows, dxy, labels)

