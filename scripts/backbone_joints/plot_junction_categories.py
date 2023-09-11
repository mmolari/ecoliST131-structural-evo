import argparse

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import pypangraph as pp
import utils as ut

from Bio import Phylo
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Given a backbone joints pangraph and the core-genome tree, it groups identical
        paths in categories and displays the categories in a figure, along with the
        tree.
        """
    )
    parser.add_argument("--pangraph", type=str)
    parser.add_argument("--tree", type=str)
    parser.add_argument("--fig", type=str)
    parser.add_argument("--len_filt", type=int, default=0)
    return parser.parse_args()


def plot_row(ax, nodes, y, colors, block_df):
    x = 0
    for node in nodes:
        color = colors[node.id]
        l = block_df["len"][node.id]
        h = block_df["duplicated"][node.id] * 0.2 + 0.4
        core = block_df["core"][node.id]
        edge_color = "black" if node.strand else "red"
        ax.barh(
            y,
            l,
            height=h,
            left=x,
            color=color,
            edgecolor=edge_color,
            hatch="." if core else None,
        )
        x += l


def plot_categories(path_categories, block_df, tree_file):
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
    N_blocks = len(block_df)
    colors = mpl.cm.get_cmap("rainbow")(np.linspace(0, 1, N_blocks))
    np.random.shuffle(colors)
    block_color = {block: colors[i] for i, block in enumerate(block_df.index)}

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

    return fig, axs


if __name__ == "__main__":
    args = parse_args()

    # extract edge name
    edge = args.pangraph.split("/")[-1].split(".")[0]

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)
    bdf = pan.to_blockstats_df()
    block_len = bdf["len"].to_dict()
    is_core = bdf["core"].to_dict()
    is_dupl = bdf["duplicated"].to_dict()

    # build paths
    paths = ut.pangraph_to_path_dict(pan)

    # optional: clean up paths
    if args.len_filt > 0:
        paths = ut.filter_paths(paths, lambda x: block_len[x] >= args.len_filt)
        bdf = bdf[bdf["len"] >= args.len_filt].copy()

    # subdivide paths into categories
    path_cat = ut.path_categories(paths)

    # plot
    fig, axs = plot_categories(path_cat, bdf, args.tree)
    axs[1].set_title(f"{edge}")
    plt.tight_layout()
    plt.savefig(args.fig, facecolor="w")
    plt.close(fig)
