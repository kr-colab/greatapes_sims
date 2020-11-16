import tskit
import pyslim
import msprime
import dendropy
import numpy as np
import pandas as pd

def add_time(ts, dt):
    '''
    This function returns a tskit.TreeSequence in which `dt`
    has been added to the times in all nodes.
    '''
    tables = ts.tables
    nodes_dict = tables.nodes.asdict()
    nodes_dict['time'] = nodes_dict['time'] + dt
    tables.nodes.set_columns(**nodes_dict)
    migrations_dict = tables.migrations.asdict()
    migrations_dict['time'] = migrations_dict['time'] + dt
    tables.migrations.set_columns(**migrations_dict)
    mutations_dict = tables.mutations.asdict()
    if not np.any(np.isnan(mutations_dict['time'])):
        mutations_dict['time'] = mutations_dict['time'] + dt
        tables.mutations.set_columns(**mutations_dict)
    return pyslim.SlimTreeSequence.load_tables(tables)


def match_nodes(tseqs, split_time):
    """
    Given two SLiM tree sequences, returns a dictionary relating
    the id in ts2 (key) to id in ts1 (item) for  node IDs in the
    two tree sequences that refer to the same node. If split time
    in ts2 (T2) is given, then only nodes before the split are
    considered. Note the only check of equivalency is the slim_id
    of the nodes.
    """
    node_mapping = np.full(tseqs[1].num_nodes, tskit.NULL)
    sids0 = np.array([n.metadata["slim_id"] for n in tseqs[0].nodes()])
    sids1 = np.array([n.metadata["slim_id"] for n in tseqs[1].nodes()])
    alive_before_split1 = tseqs[1].tables.nodes.time >= split_time
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

def sub_metadata(tseqs):
    """
    Work around current bug in `tskit.union`: subbing top-level metadata
    so they match.
    """
    tables0 = tseqs[0].tables
    tables0.metadata = tseqs[1].tables.metadata
    tseqs[0] = pyslim.SlimTreeSequence.load_tables(tables0)

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
    return msprime.RateMap(breaks, rates)


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
        assert subset.shape[0] == 1
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
    for node in tree.postorder_node_iter(filter_fn = lambda node: node.is_internal()):
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
                tseqs.append(pyslim.load(trees_path+child.taxon.label+"_"+rand_id+"_rep"+rep+".trees"))
            else:
                tseq, p = in_tseqs.pop(child.taxon.label)
                tseqs.append(tseq)
                pops += p
                del tseq
        #check if times need be shifted
        print(f"Before shift\ttime 0: {tseqs[0].max_root_time}\ttime 1: {tseqs[1].max_root_time}")
        if history_len[1] > history_len[0]:
            tseqs[0] = add_time(tseqs[0], history_len[1]-history_len[0])
        elif history_len[0] > history_len[1]:
            tseqs[1] = add_time(tseqs[1], history_len[0]-history_len[1])
        print(f"After shift\ttime 0: {tseqs[0].max_root_time}\ttime 1: {tseqs[1].max_root_time}")
        node_mapping = match_nodes(tseqs, node.age)
        sub_metadata(tseqs)
        in_tseqs[node.taxon.label] = (tseqs[0].union(tseqs[1], node_mapping), pops)
    assert len(in_tseqs) == 1
    return in_tseqs[list(in_tseqs.keys())[0]]
