# %%

import json
import pathlib

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import utils as ut

import pypangraph as pp

from Bio import Phylo
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

    # sort by count
    path_cat = [(count, nodes[path], iso_list[path]) for path, count in n_paths.items()]
    path_cat.sort(key=lambda x: x[0], reverse=True)
    return path_cat


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

# # plot paths

# tree = Phylo.read("data/named_tree.nwk", "newick")
# tree.ladderize()
# str_order = [n.name for n in tree.get_terminals()]


# def plot_row(ax, nodes, i, colors, block_df):
#     x = 0
#     for node in nodes:
#         color = colors[node.id]
#         l = block_df["len"][node.id]
#         w = block_df["dupl"][node.id] * 3 + 3
#         ax.plot([x, x + l], [i, i], color=color, linewidth=w)
#         x += l


# def plot_gap(ax, Js, e, str_order, lengths):
#     cmap = mpl.cm.get_cmap("rainbow")
#     colors = defaultdict(lambda: cmap(np.random.rand()))

#     for i, iso in enumerate(str_order):
#         j = Js[iso]
#         if j.left != e.left:
#             j = j.invert()
#             assert j.left == e.left
#         plot_row(ax, j.center.nodes, i, colors, lengths)

#     ax.set_yticks(np.arange(len(str_order)))
#     ax.set_yticklabels(str_order)
#     ax.set_xlabel("position")


# svfld = pathlib.Path("figs") / "core_junctions"
# svfld.mkdir(exist_ok=True)

# for e, j_dict in junctions.items():
#     nn_j = len([j for j in j_dict.values() if len(j.center.nodes) > 0])
#     if nn_j < 2:
#         continue
#     fig, ax = plt.subplots(figsize=(20, 20))
#     plot_gap(ax, j_dict, e, str_order, block_cts[e])
#     plt.title(f"{e}")
#     plt.tight_layout()
#     plt.savefig(
#         svfld / f"{e.left.to_str_id()}__{e.right.to_str_id()}.png", facecolor="w"
#     )
#     plt.show()

# %%


def plot_row(ax, nodes, y, colors, block_df):
    x = 0
    for node in nodes:
        color = colors[node.id]
        l = block_df["len"][node.id]
        h = block_df["dupl"][node.id] * 0.2 + 0.4
        edge_color = "black" if node.strand else "red"
        ax.barh(
            y,
            l,
            height=h,
            left=x,
            color=color,
            edgecolor=edge_color,
            # lw=3,
        )
        x += l


def plot_categories(e, path_categories, block_df, tree_file):
    # load tree
    tree = Phylo.read(tree_file, "newick")
    tree.ladderize()
    leaves = [n for n in tree.get_terminals()]

    # assign colors to leaves
    C = len(path_categories)

    if C <= 10:
        path_colors = mpl.cm.get_cmap("tab10")(np.arange(C))
    elif C <= 20:
        path_colors = mpl.cm.get_cmap("tab20")(np.arange(C))
    else:
        path_colors = mpl.cm.get_cmap("jet")(np.linspace(0, 1, C))

    strain_color = defaultdict(lambda: "white")
    for i, (_, _, isolates) in enumerate(path_categories):
        for iso in isolates:
            strain_color[iso] = path_colors[i]

    # assign color to blocks
    cmap = mpl.cm.get_cmap("rainbow")
    block_color = defaultdict(lambda: cmap(np.random.rand()))

    fig, axs = plt.subplots(
        1, 2, figsize=(20, 15), gridspec_kw={"width_ratios": [1, 3]}
    )

    ax = axs[0]
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda x: x.name if x in leaves else "",
        label_colors=lambda x: strain_color[x],
    )
    ax.set_title(f"{e}")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax = axs[1]
    for i, (count, nodes, isolates) in enumerate(path_categories):
        plot_row(ax, nodes, -i, block_color, block_df)
        ax.text(
            0,
            -i + 0.45,
            f"path {i+1} | n = {count}",
            color=path_colors[i],
            fontsize=20,
        )

    ax.set_yticks([])
    ax.set_ylim(-len(path_categories) + 0.5, 0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("position")

    plt.tight_layout()
    return fig, ax


svfld = pathlib.Path("figs") / "core_junctions_categories"
svfld.mkdir(exist_ok=True)

for e, pc in path_cts.items():
    if len(pc) == 0:
        continue
    if (len(pc) == 1) and (pc[0][0] < 2):
        continue
    fig, ax = plot_categories(e, pc, block_cts[e], "data/named_tree.nwk")
    plt.savefig(
        svfld / f"{e.left.to_str_id()}__{e.right.to_str_id()}.png", facecolor="w"
    )
    plt.close(fig)

    # break
# %%
