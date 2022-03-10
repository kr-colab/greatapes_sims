import tskit
import pyslim
import msprime
import dendropy
import numpy as np
import pandas as pd
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

def get_slim_id(node, mschema=pyslim.slim_metadata_schemas["node"]):
    # gets slim id for a node
    # if treeseq didnt have metadata schemas, uses pyslim
    sid = None
    try:
        sid = node.metadata['slim_id']
    except:
        sid = mschema.decode_row(node.metadata)['slim_id']
    return sid

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
    sids0 = np.array([get_slim_id(node) for node in ts1.nodes()])
    sids1 = np.full(ts2.num_nodes, -1)
    alive_before_split1 = np.full(ts2.num_nodes, False)
    for i, node in enumerate(ts2.nodes()):
        sids1[i] = get_slim_id(node)
        alive_before_split1[i] = (node.time >= split_time)
    assert np.all(sids1!=-1)
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


def subtree(focal, edges, taxon_namespace, nodes = None, child="edge", parent="parent"):
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
        print(node.taxon.label)
        subset = meta[(meta.edge==node.taxon.label) & (meta.rand_id == rand_id)]
        print(subset)
        subset.drop(columns=["date"], inplace=True)
        assert len(subset.drop_duplicates()) == 1
        n_gens = np.floor(subset.gens.values[0]/subset.rescf.values[0])
        node.edge_length= n_gens
        #print(node.edge_length)
        #print(node.distance_from_tip())
    tree.calc_node_root_distances(return_leaf_distances_only=False)
    tree.calc_node_ages(ultrametricity_precision=False, is_force_max_age=True)
    return tree

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
        print("Garbage collector: collected", "%d objects." % collected, flush=True)
        assert len(node.child_nodes()) == 2, "Polytomies are not supported."
        tseqs = []
        pops = []
        history_len = []
        print(node.taxon.label, "\t", node.age, sep="", flush=True)
        for child in node.child_nodes():
            print("\t"+child.taxon.label+"\t"+str(child.root_distance)+"\t"+str(child.age), flush=True)
            history_len.append(child.root_distance+child.age)
            if child.is_leaf():
                pops.append(child.taxon.label)
                print(trees_path+child.taxon.label+"_"+rand_id+"_rep"+rep+".trees", flush=True)
                tpath = trees_path+child.taxon.label+"_"+rand_id+"_rep"+rep+".trees"
                tseqs.append(tskit.load(tpath))
                print(gc.get_stats(), flush=True)
                collected = gc.collect()
            else:
                tseq, p = in_tseqs.pop(child.taxon.label)
                tseqs.append(tseq)
                pops += p
                del tseq
        assert len(tseqs) == 2
        #check if times need be shifted
        print(f"Before shift\ttime 0: {tseqs[0].max_root_time}\ttime 1: {tseqs[1].max_root_time}", flush=True)
        if history_len[1] > history_len[0]:
            tseqs[0] = refactor_time(tseqs[0], history_len[1]-history_len[0])
        elif history_len[0] > history_len[1]:
            tseqs[1] = refactor_time(tseqs[1], history_len[0]-history_len[1])
        print(f"After shift\ttime 0: {tseqs[0].max_root_time}\ttime 1: {tseqs[1].max_root_time}", flush=True)
        print("Matching nodes", flush=True)
        node_mapping = match_nodes(tseqs[0], tseqs[1], node.age)
        print("Union\'ing pops: ", pops, flush=True)
        tsu = tseqs[0].union(tseqs[1], node_mapping)
        in_tseqs[node.taxon.label] = tsu, pops
    assert len(in_tseqs) == 1
    assert node.taxon.label == list(in_tseqs.keys())[0]
    return (tsu, pops)

def sample_from_ts(ts, sample_size=None, contemporary=True, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    # Returns a dicionary with `sample_size` samples for each population.
    # If `contemporary` is True, only samples whose times are close to
    # the minimum times in the sample's population are kept.
    # Useful if tree has multiple pops and is not ultrametric.
    min_times = {pop.id: np.inf for pop in ts.populations()}
    if contemporary:
        for node in ts.nodes():
            if node.flags != tskit.NODE_IS_SAMPLE:
                continue
            if node.time < min_times[node.population]:
                min_times[node.population] = node.time
    samples = {}
    for node in ts.nodes():
        if node.flags != tskit.NODE_IS_SAMPLE:
            continue
        if node.population not in samples:
            samples[node.population] = []
        if np.isclose(node.time, min_times[node.population]) or (not contemporary):
            samples[node.population].append(node.id)
    if sample_size is not None:
        for pop, sample in samples.items():
            samples[pop] = rng.choice(sample, sample_size, replace=False)
    return samples
