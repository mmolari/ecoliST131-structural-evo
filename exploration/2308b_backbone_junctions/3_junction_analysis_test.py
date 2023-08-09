# %%

import json

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import utils as ut

import pypangraph as pp

from collections import Counter, defaultdict

LEN_THR = 200


def load_junctions(j_file):
    # load json file
    with open(j_file, "r") as f:
        str_junctions = json.load(f)

    # from strings to objects
    junctions = {}
    for e_str, d in str_junctions.items():
        e = ut.Edge.from_str_id(e_str)
        junctions[e] = {}
        for iso, j_str in d.items():
            junctions[e][iso] = ut.Junction.from_list(j_str)

    return junctions


# %%

# load pangraph
pan = pp.Pangraph.load_json("data/pangraph.json")
bdf = pan.to_blockstats_df()
block_len = bdf["len"].to_dict()

# %%

# load junctions
junctions = load_junctions("data/junctions.json")


def keep_node(node):
    keep = block_len[node.id] >= LEN_THR
    return keep


def filter_paths(j_dict, keep_f):
    """Filter out nodes using the input filter function."""
    for iso, j in j_dict.items():
        j.center.nodes = list(filter(keep_f, j.center.nodes))


def path_block_count_df(j_dict):
    """Count occurrences of blocks"""
    count = Counter()
    n_strains = Counter()

    for iso, j in j_dict.items():
        bids = [n.id for n in j.center.nodes]
        count.update(bids)
        n_strains.update(set(bids))

    df = pd.DataFrame.from_dict(count, orient="index", columns=["count"])
    df["n_strains"] = df.index.map(n_strains)
    df["dupl"] = df["count"] > df["n_strains"]

    return df


def path_edge_count(j_dict):
    """Count internal edges of paths"""
    edge_ct = Counter()
    for iso, j in j_dict.items():
        for n1, n2 in zip(j.center.nodes[:-1], j.center.nodes[1:]):
            edge_ct.update([ut.Edge(n1, n2)])
    return dict(edge_ct)


def find_mergers(edge_ct, block_ct):
    """Create a dictionary source -> sinks of block-ids to be merged"""
    mergers = {}
    for e, ec in edge_ct.items():
        bl, br = e.left.id, e.right.id
        if (ec == block_ct[bl]) and (ec == block_ct[br]):
            # merge
            if bl in mergers:
                mergers[br] = mergers[bl]
            elif br in mergers:
                mergers[bl] = mergers[br]
            else:
                mergers[br] = bl
    return mergers


def perform_mergers(mergers, j_dict, block_ct):
    """Performs selected mergers in paths and block stats"""
    for source, sink in mergers.items():
        # modify block count dataframe
        block_ct.loc[sink, "len"] += block_ct["len"][source]
        assert block_ct["count"][sink] == block_ct["count"][source]
        assert block_ct["n_strains"][sink] == block_ct["n_strains"][source]
        block_ct.drop(source, inplace=True)

    def keep_node(node):
        return node.id not in mergers.keys()

    filter_paths(j_dict, keep_node)


def path_categories(anchor_edge, j_dict):
    """Returns a list of touples, one per non-empty path, with the following info:
    (count, path, [list of isolates])"""
    iso_list = defaultdict(list)
    n_paths = defaultdict(int)
    nodes = {}
    for iso, j in j_dict.items():
        if j.left != anchor_edge.left:
            j = j.invert()
        assert j.left == anchor_edge.left
        assert j.right == anchor_edge.right

        path = j.center
        if len(path.nodes) > 0:
            n_paths[path] += 1
            iso_list[path].append(iso)
            nodes[path] = path.nodes

    return [(count, nodes[path], iso_list[path]) for path, count in n_paths.items()]


block_cts = {}
path_cts = {}
for e, j_dict in junctions.items():
    filter_paths(j_dict, keep_node)
    # check length threshold
    assert np.all(
        [
            np.all([block_len[n.id] >= LEN_THR for n in j.center.nodes])
            for j in j_dict.values()
        ]
    )

    # cound number of blocks
    block_ct = path_block_count_df(j_dict)
    block_ct["len"] = block_ct.index.map(block_len)

    # count number of edges
    edge_ct = path_edge_count(j_dict)

    # merge transitive edges
    mergers = find_mergers(edge_ct, block_ct["count"].to_dict())
    perform_mergers(mergers, j_dict, block_ct)

    # TODO: paralog split?

    block_cts[e] = block_ct
    path_cts[e] = path_categories(e, j_dict)

# %%

# plot paths

import pathlib
from Bio import Phylo
from collections import defaultdict

tree = Phylo.read("data/named_tree.nwk", "newick")
tree.ladderize()
str_order = [n.name for n in tree.get_terminals()]


def plot_row(ax, nodes, i, colors, block_df):
    x = 0
    for node in nodes:
        color = colors[node.id]
        l = block_df["len"][node.id]
        w = block_df["dupl"][node.id] * 3 + 3
        ax.plot([x, x + l], [i, i], color=color, linewidth=w)
        x += l


def plot_gap(ax, Js, e, str_order, lengths):
    cmap = mpl.cm.get_cmap("rainbow")
    colors = defaultdict(lambda: cmap(np.random.rand()))

    for i, iso in enumerate(str_order):
        j = Js[iso]
        if j.left != e.left:
            j = j.invert()
            assert j.left == e.left
        plot_row(ax, j.center.nodes, i, colors, lengths)

    ax.set_yticks(np.arange(len(str_order)))
    ax.set_yticklabels(str_order)
    ax.set_xlabel("position")


svfld = pathlib.Path("figs") / "core_junctions"
svfld.mkdir(exist_ok=True)

for e, j_dict in junctions.items():
    nn_j = len([j for j in j_dict.values() if len(j.center.nodes) > 0])
    if nn_j < 2:
        continue
    fig, ax = plt.subplots(figsize=(20, 20))
    plot_gap(ax, j_dict, e, str_order, block_cts[e])
    plt.title(f"{e}")
    plt.tight_layout()
    plt.savefig(
        svfld / f"{e.left.to_str_id()}__{e.right.to_str_id()}.png", facecolor="w"
    )
    plt.show()

# %%
