# %%
import numpy as np
import pandas as pd
import pypangraph as pp
from Bio import Phylo
import matplotlib as mpl
import matplotlib.pyplot as plt
import pathlib
from collections import defaultdict
import argparse

left_clr = "#2db7a7"
right_clr = "#cf4b36"


def parse_args():
    parser = argparse.ArgumentParser(description="Draw linear junctions")
    parser.add_argument("--hs", type=str, help="Hotspot name")
    return parser.parse_args()


def show():
    # plt.show()
    plt.close()


def despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def load_data(hs):
    fname = f"data/{hs}/joint_graph.json"
    pan = pp.Pangraph.load_json(fname)

    fname = "data/tree.nwk"
    tree = Phylo.read(fname, "newick")

    MGEs = pd.read_csv(f"res/{hs}/tool_annotations.csv")

    hh_genes = pd.read_csv("data/hh_genes.csv")
    hh_genes = hh_genes[
        ["Gene symbol of upstream gene", "Gene symbol of downstream gene"]
    ]
    hh_genes.columns = ["up", "down"]

    ann = pd.read_csv(f"res/{hs}/gbk_annotations.csv")

    return pan, tree, MGEs, hh_genes, ann


def color_generator():
    cm = plt.get_cmap("tab20")
    for i in range(20):
        yield cm(i)
    cm = plt.get_cmap("tab20b")
    for i in range(20):
        yield cm(i)
    cm = plt.get_cmap("tab20c")
    for i in range(20):
        yield cm(i)
    cm = plt.get_cmap("hsv")
    while True:
        yield cm(np.random.rand())


def bandage_color_df(block_colors):
    colour_df = pd.DataFrame.from_dict(block_colors, orient="index", columns=["Colour"])
    colour_df["Colour"] = colour_df["Colour"].apply(
        lambda x: mpl.colors.to_hex(x) if not isinstance(x, str) else x,
    )
    colour_df.index.name = "Name"
    return colour_df


def figure_setup(tree):
    fig, axs = plt.subplots(
        1,
        2,
        figsize=(11.69, 8.27),
        sharey=True,
        gridspec_kw={"width_ratios": [1, 10]},
    )
    ax = axs[0]
    Phylo.draw(
        tree, axes=ax, label_func=lambda x: "", do_show=False, show_confidence=False
    )
    return fig, axs


def display_MGEs(ax, MGEs, tree):
    strain_y = {l.name: k + 1 for k, l in enumerate(tree.get_terminals())}

    colors = {
        "ISEScan": "C0",
        "defensefinder": "C2",
        "integronfinder": "C3",
        "genomad": "C4",
    }

    for _, row in MGEs.iterrows():
        strain = row["iso"]
        y = strain_y[strain]
        s, e = row["start"], row["end"]
        kind = row["kind"]
        c = colors[kind]
        ax.plot([s, e], [y, y], color=c, lw=1)


def linear_junction_draw(tree, pan, MGEs=None, color_blocks=True):

    fig, axs = figure_setup(tree)

    strains = [cl.name for cl in tree.get_terminals()]

    if color_blocks:
        gen = color_generator()
        block_colors = defaultdict(lambda: next(gen))
    else:
        block_colors = defaultdict(lambda: "silver")

    # assign colors to first/last blocks
    Bs = pan.paths[strains[0]].block_ids
    block_colors[Bs[0]] = left_clr
    block_colors[Bs[-1]] = right_clr

    ax = axs[1]
    for k, iso in enumerate(strains):
        y = k + 1
        if iso not in pan.strains():
            continue
        path = pan.paths[iso]
        Bs = path.block_ids
        Ss = path.block_strands
        Ps = path.block_positions
        for n_block in range(len(Bs)):
            bid = Bs[n_block]
            strand = Ss[n_block]
            x0, x1 = Ps[n_block], Ps[n_block + 1]
            color = block_colors[bid]
            edgecolor = "none"
            if color_blocks and not strand:
                edgecolor = "k"
            ax.barh(
                y,
                x1 - x0,
                left=x0,
                height=0.8,
                color=color,
                edgecolor=edgecolor,
                linewidth=0.5,
            )
    ax.set_xlabel("base pairs")

    # optional MGE display
    if MGEs is not None:
        display_MGEs(ax, MGEs, tree)

    # despine
    for ax in axs:
        despine(ax)

    plt.tight_layout()

    block_colors = bandage_color_df(block_colors)

    return fig, block_colors


def add_hh_genes(ax, tree, ann, hh_genes):

    up = set(hh_genes["up"])
    down = set(hh_genes["down"])

    # exclude NaNs:
    up = {g for g in up if isinstance(g, str)}
    down = {g for g in down if isinstance(g, str)}

    strain_y = {l.name: k + 1 for k, l in enumerate(tree.get_terminals())}
    mask = ann["gene"].isin(up | down)
    mask &= ann["kind"] == "CDS"

    for _, row in ann[mask].iterrows():
        b, e = row["start"], row["end"]
        strand = row["strand"]
        y = strain_y[row["iso"]]
        gene = row["gene"]
        if gene in up:
            clr = "blue"
        else:
            clr = "red"
        marker = ">" if strand == 1 else "<"
        ax.scatter((b + e) / 2, y, color=clr, marker=marker, s=10, zorder=10)

    # add legend
    genes = ann[mask]["gene"].unique()
    for g in genes:
        if g in up:
            clr = "blue"
        else:
            clr = "red"
        ax.scatter([], [], color=clr, label=g, marker=">", s=10)
    ax.legend(loc="upper right", title="genes", title_fontsize=10, fontsize=8)


if __name__ == "__main__":

    args = parse_args()
    hs = args.hs
    # hs = "CIRMBUYJFK_f__CWCCKOQCWZ_r"

    pan, tree, MGEs, hh_genes, ann = load_data(hs)

    fig_fld = pathlib.Path(f"figs/{hs}")
    fig_fld.mkdir(parents=True, exist_ok=True)

    res_fld = pathlib.Path(f"res/{hs}")
    res_fld.mkdir(parents=True, exist_ok=True)

    bdf = pan.to_blockstats_df()

    # fig, block_colors = linear_junction_draw(tree, pan)
    # plt.savefig(fig_fld / f"linear_block_repr.png", dpi=300)
    # plt.savefig(fig_fld / f"linear_block_repr.svg")
    # plt.show()

    # block_colors.to_csv(res_fld / "export/block_colours.csv")

    fig, block_colors = linear_junction_draw(tree, pan, MGEs=MGEs, color_blocks=False)
    add_hh_genes(fig.get_axes()[1], tree, ann, hh_genes)
    plt.savefig(fig_fld / f"linear_block_repr_MGEs.png", dpi=300)
    plt.savefig(fig_fld / f"linear_block_repr_MGEs.svg")
    plt.show()

# %%
