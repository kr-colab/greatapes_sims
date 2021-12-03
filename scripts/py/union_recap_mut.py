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
#from helper_functions import *
import operator
import gc

def refactor_time(ts, factor, operation=operator.iadd):
    '''
    This function returns a tskit.TreeSequence in which the times columns in
    the nodes, migrations and mutations table have been refactored by `factor`
    using an `operation`.
    `operation` is a function from the module operator (e.g., `iadd` or `imul`)
    `operator.iadd(x, factor)` is equivalent to `x += factor`
    `operator.imul(x, factor)` is equivalent to `x *= factor`
    '''
    tables = ts.dump_tables()
    for table in (tables.migrations, tables.mutations, tables.nodes):
        if not np.any(np.isnan(table.time)):
            table.time = operation(table.time, factor)
    return tables.tree_sequence()

def match_nodes(ts1, ts2, split_time):
    """
    Given two SLiM tree sequences, returns a dictionary relating
    the id in ts2 (key) to id in ts1 (item) for  node IDs in the
    two tree sequences that refer to the same node. If split time
    in ts2 (T2) is given, then only nodes before the split are
    considered. Note the only check of equivalency is the slim_id
    of the nodes.
    """
    node_mapping = np.full(ts2.num_nodes, tskit.NULL)
    sids0 = np.array([n.metadata["slim_id"] for n in ts1.nodes()])
    sids1 = np.array([n.metadata["slim_id"] for n in ts2.nodes()])
    alive_before_split1 = ts2.tables.nodes.time >= split_time
    sorted_ids0 = np.argsort(sids0)
    matches = np.searchsorted(
        sids0,
        sids1,
        side='left',
        sorter=sorted_ids0)
    is_1in0 = np.isin(sids1, sids0)
    both = np.logical_and(alive_before_split1, is_1in0)
    node_mapping[both] = sorted_ids0[matches[both]]
    return node_mapping

def msp_mutation_rate_map(intervals, total_rate, intervals_rate, length):
    """
    Takes a `pd.DataFrame` with three columns (?, start, end), with 0-indexed [start, end) intervals.
    Returns breaks and rates to use with `msprime.mutate`, in which the rate for `msprime` will be
    `total_rate-intervals_rate` within the intervals.
    """
    breaks = [0]
    rates = []
    for (i, c, start, end) in intervals.itertuples():
        if start not in breaks:
            breaks.append(start)
            rates.append(total_rate)
        breaks.append(end)
        rates.append(total_rate-intervals_rate)
    if not np.isclose(breaks[-1], length):
        breaks.append(length)
        rates.append(total_rate)
    return msprime.RateMap(position=breaks, rate=rates)


def subtree(focal, edges, taxon_namespace, nodes = None):
    """
    Returns a dictionary of `dendropy.Node` objects from a `pandas.DataFrame`
    with two columns: `edge` and `parent`, which specifies
    the edge-parent relationships. Only nodes below focal are returned.
    """
    if nodes == None:
        nodes = {}
    if not focal in nodes:
        nodes[focal] = dendropy.Node(taxon=taxon_namespace.get_taxon(focal))
    for i, row in edges.iterrows():
        if row.parent == focal:
            if not row.edge in nodes:
                nodes[row.edge] = dendropy.Node(taxon=taxon_namespace.get_taxon(row.edge))
            nodes[focal].add_child(nodes[row.edge])
            nodes = subtree(row.edge, edges, taxon_namespace, nodes)
    return nodes

def build_tree_from_df(edges):
    """
    Returns a `dendropy.Tree` from a `pandas.DataFrame` with edge-parent
    relationships.
    """
    root_name = edges.edge[edges.parent==""][0]
    taxon_namespace = dendropy.TaxonNamespace(edges.edge.values.tolist())
    nodes = subtree(root_name, edges, taxon_namespace)
    tree = dendropy.Tree(seed_node = nodes[root_name], taxon_namespace=taxon_namespace)
    tree.ladderize(ascending=False)
    return(tree)

def add_blen_from_meta(tree, meta, rand_id):
    """
    `meta` is a `pandas.DataFrame` with columns `edge`, `rand_id`, `gens` and
    `rescf`. This function adds branch lengths to the `dendropy.Tree` object
    using the info in the `meta`.
    """
    rep = '0' # all reps should be the have the sam blens anyway!
    # traversing through tree -- annotating lengths
    for node in tree.postorder_node_iter():
        #print(node.taxon.label)
        subset = meta[(meta.edge==node.taxon.label) & (meta.rand_id == rand_id)]
        #print(subset)
        subset.drop(columns=["date"], inplace=True)
        assert len(subset.drop_duplicates()) == 1
        n_gens = np.floor(subset.gens.values[0]/subset.rescf.values[0])
        node.edge_length= n_gens
        #print(node.edge_length)
        #print(node.distance_from_tip())
    tree.calc_node_root_distances(return_leaf_distances_only=False)
    tree.calc_node_ages(ultrametricity_precision=False, is_force_max_age=True)
    return tree

@profile
def pyslimload(x):
    return pyslim.load(x)

@profile
def union_tseqs(tree, rand_id, rep, trees_path):
    """
    Given a `dendropy.tree` object with annotated `edge_lengths`, a `rand_id`
    identifier and a replicate number `rep`, this performs the
    `tskit.TableCollection.union` of all leaves in the phylogenetic tree.
    """
    in_tseqs = {}
    tsu = None
    for node in tree.postorder_node_iter(filter_fn = lambda node: node.is_internal()):
        del tsu
        collected = gc.collect()
        print("Garbage collector: collected", "%d objects." % collected)
        assert len(node.child_nodes()) == 2, "Polytomies are not supported."
        tseqs = []
        pops = []
        history_len = []
        print(node.taxon.label, "\t", node.age, sep="")
        for child in node.child_nodes():
            print("\t"+child.taxon.label+"\t"+str(child.root_distance)+"\t"+str(child.age))
            history_len.append(child.root_distance+child.age)
            if child.is_leaf():
                pops.append(child.taxon.label)
                print(trees_path+child.taxon.label+"_"+rand_id+"_rep"+rep+".trees")
                tseqs.append(trees_path+child.taxon.label+"_"+rand_id+"_rep"+rep+".trees")
            else:
                tseq, p = in_tseqs.pop(child.taxon.label)
                tseqs.append(tseq)
                pops += p
                del tseq
        assert len(tseqs) == 2
        #check if times need be shifted
        collected = gc.collect()
        print("Before pyslim.load garbage collector: collected", "%d objects." % collected)
        ts1 = pyslimload(tseqs[0])
        collected = gc.collect()
        print("After pyslim.load garbage collector: collected", "%d objects." % collected)
        ts2 = pyslimload(tseqs[1])
        collected = gc.collect()
        print("Garbage collector: collected", "%d objects." % collected)
        print(f"Before shift\ttime 0: {ts1.max_root_time}\ttime 1: {ts2.max_root_time}")
        if history_len[1] > history_len[0]:
            ts1 = refactor_time(ts1, history_len[1]-history_len[0])
        elif history_len[0] > history_len[1]:
            ts2 = refactor_time(ts2, history_len[0]-history_len[1])
        print(f"After shift\ttime 0: {ts1.max_root_time}\ttime 1: {ts2.max_root_time}")
        print("Matching nodes")
        node_mapping = match_nodes(ts1, ts2, node.age)
        #sub_metadata(tseqs)
        print("Union\'ing pops: ", pops)
        tsu = ts1.union(ts2, node_mapping)
        del ts1
        del ts2
        tsuname = trees_path+node.taxon.label+"_"+rand_id+"_rep"+rep+".trees"
        tsu.dump(tsuname)
        in_tseqs[node.taxon.label] = tsuname, pops
    assert len(in_tseqs) == 1
    assert node.taxon.label == list(in_tseqs.keys())[0]
    return (tsu, pops)

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
sims_full.columns = header.columns
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
print("Loading union-ing tree sequences!")
tsu,  pops = union_tseqs(tree,args["rand_id"],args["rep"], trees_path+args['rand_id']+"/")
tcu = tsu.dump_tables()
del tsu
print("Fixing metadata")
if np.any(np.isnan(tcu.mutations.time)):
    # TODO: remove this once using slim that adds time
    print("SLiM tree sequence did not have mutation times.")
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

print("Dumping union\'ed treeseq")
tsu.dump(union_path)
with open(pops_path, "w") as f:
    f.write(str(pops))

# recapitating
# doing a hacky recapitation bc pyslim hasn't been updated to use msprime1.0
print("Recapitating")
demography = msprime.Demography.from_tree_sequence(tsu)
for pop in demography.populations:
    pop.initial_size=args["recapN"]
recomb_map = msprime.RateMap.read_hapmap(rec_hap_path, position_col=1, rate_col=2)
recap_tsu = msprime.sim_ancestry(initial_state=tsu, recombination_rate=recomb_map, demography=demography)
del tsu # too much ram
print(slim_gen, recap_tsu.max_root_time, recap_tsu.num_mutations)

# mutating
mut_map = msp_mutation_rate_map(exons, args["total_mut_rate"], args["region_mut_rate"], int(recap_tsu.sequence_length))
# figuring out the max id slim used
print("Mutating!")
max_id = -1
for mut in recap_tsu.mutations():
    for d in mut.derived_state.split(","):
        max_id = max(max_id, int(d))
model = msprime.SLiMMutationModel(type=3, next_id=max_id+1)
print("Before mutate:", recap_tsu.num_mutations)
# adding mutations to recapitated part
recap_tsu = msprime.sim_mutations(recap_tsu, start_time=slim_gen, model=model, rate=args["total_mut_rate"], keep=True)
print("Mutations added in the recapitation:", recap_tsu.num_mutations)
# adding mutations to the SLiM part
recap_tsu = msprime.sim_mutations(recap_tsu, end_time=slim_gen, model=model, rate=mut_map, keep=True)
print("Total mutations:", recap_tsu.num_mutations)
recap_tsu.dump(recap_mut_path)
