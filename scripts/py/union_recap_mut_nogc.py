#!/usr/bin/env python
import tskit
import msprime
import dendropy
import glob
import numpy as np
import pandas as pd
import warnings
import functools
import argparse
import operator
import os
import pyslim
from helper_functions import *
import operator

# variables
parser = argparse.ArgumentParser(description='Perform the union and output stats from tree sequences.')
parser.add_argument('rand_id', type=str)
parser.add_argument('rep', type=str)
parser.add_argument('total_mut_rate', type=float)
parser.add_argument('region_mut_rate', type=float)
parser.add_argument('--rescf', type=int, default=1)
parser.add_argument('--recapN', type=int, default=10000, required=False)
parser.add_argument('--seed', type=int, default=8991, required=False)

# 7TMN8S0F7TVAI2C 0 2e-08 1.2001e-08 --rescf 1 --recapN 10000 --seed 943284230
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
rec_hap_path = f"../../output/maps/{args['rand_id']}_recrate.hapmap"
ex_path = f"../../output/maps/{args['rand_id']}_exons.tsv"


edges_path = "../../data/meta/edges_meta.tsv"
sims_full_path = "../../output/sims_info.tsv"
sims_header_path = "../../output/header_sims_info.tsv"
trees_path = "../../output/"

## loading metadata
# edges contains all the edges and info about N and number of generations
print("Loading metadata stuff")
edges = pd.read_csv(edges_path,sep="\t")
edges.parent= edges.parent.fillna("")
edges["edge"] = edges["edge"].str.replace('_','-')
edges["parent"] = edges["parent"].str.replace('_','-')
# sims_full relate args["rand_id"]s to simulation parameters, for instance
#generations per branch in the phylo treee
sims_full= pd.read_csv(sims_full_path,sep="\t", header=None)
header = pd.read_csv(sims_header_path,sep="\t")
cols = header.columns.to_list()
if cols[-1] != "poscoefdecayeps":
    cols.append("poscoefdecayeps")
sims_full.columns = cols
# exons file
exons = pd.read_csv(ex_path,sep="\t")
# removing extraneous columns?
exons = exons.iloc[:,:3]

# getting all output files and grouping by args["rand_id"] and args["rep"]
tree_files = glob.glob(trees_path+args['rand_id']+"/*[0-9].trees")
pattern = f"{args['rand_id']}_rep{args['rep']}"
n_matches = sum(1 for file in tree_files if pattern in file)
# making sure we got all the files
assert n_matches == edges.shape[0]


# getting the phylo tree adn annotating with branch lengths,
tree = build_tree_from_df(edges)
tree = add_blen_from_meta(tree, sims_full, args["rand_id"])

# performing the union
union_path = f"{trees_path}{args['rand_id']}/{args['rand_id']}_rep{args['rep']}.union.trees"
recap_mut_path = f"{trees_path}{args['rand_id']}/{args['rand_id']}_rep{args['rep']}.union.recap.mut.trees"
pops_path = f"{trees_path}{args['rand_id']}/{args['rand_id']}_rep{args['rep']}.pops"
print("Loading union-ing tree sequences!", flush=True)
tsu,  pops = union_tseqs(tree,args["rand_id"],args["rep"], trees_path+args['rand_id']+"/")
print(pops, flush=True)
tcu = tsu.dump_tables()
del tsu
if np.any(np.isnan(tcu.mutations.time)):
    # TODO: remove this once using slim that adds time
    print("SLiM tree sequence did not have mutation times.", flush=True)
    tcu.compute_mutation_times()
tsu = tcu.tree_sequence()
del tcu

# asserting within population coalescen
assert len(set([tsu.node(u).population for t in tsu.trees() for u in t.roots])) == 1
# adjusting times
assert tsu.max_root_time.is_integer()
slim_gen = int(tsu.max_root_time) * args['rescf']
# refactoring time if simulation was run with rescaling
if args['rescf'] > 1:
    tsu = refactor_time(tsu, args['rescf'], operator.imul)
print("Dumping union\'ed treeseq", flush=True)
tsu.dump(union_path)
with open(pops_path, "w") as f:
    f.write(str(pops))
# recapitating
# doing a hacky recapitation bc pyslim hasn't been updated to use msprime1.0
print("Recapitating", flush=True)
demography = msprime.Demography.from_tree_sequence(tsu)
for pop in demography.populations:
    pop.initial_size=args["recapN"]
recomb_map = msprime.RateMap.read_hapmap(rec_hap_path, position_col=1, rate_col=2)
recap_tsu = msprime.sim_ancestry(initial_state=tsu, recombination_rate=recomb_map, demography=demography)
del tsu # too much ram
print(slim_gen, recap_tsu.max_root_time, recap_tsu.num_mutations, flush=True)

# mutating
mut_map = msp_mutation_rate_map(exons, args["total_mut_rate"], args["region_mut_rate"], int(recap_tsu.sequence_length))
# figuring out the max id slim used
print("Mutating!", flush=True)
max_id = -1
for mut in recap_tsu.mutations():
    for d in mut.derived_state.split(","):
        max_id = max(max_id, int(d))
model = msprime.SLiMMutationModel(type=3, next_id=max_id+1)
print("Before mutate:", recap_tsu.num_mutations, flush=True)
# adding mutations to recapitated part
recap_tsu = msprime.sim_mutations(recap_tsu, start_time=slim_gen, model=model, rate=args["total_mut_rate"], keep=True)
print("Mutations added in the recapitation:", recap_tsu.num_mutations, flush=True)
# adding mutations to the SLiM part
recap_tsu = msprime.sim_mutations(recap_tsu, end_time=slim_gen, model=model, rate=mut_map, keep=True)
print("Total mutations:", recap_tsu.num_mutations, flush=True)
recap_tsu.dump(recap_mut_path)
