#  %%
# draw paths with gene annotations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
import pypangraph as pp
import argparse
import pathlib
from collections import defaultdict

left_clr = "#2db7a7"
right_clr = "#cf4b36"


def parse_args():
    parser = argparse.ArgumentParser(description="Draw linear junctions")
    parser.add_argument("--hs", type=str, help="junction name")
    return parser.parse_args()


def show():
    # plt.show()
    plt.close()


def despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def load_data(hs):
    # load tree
    tree = Phylo.read(f"data/tree.nwk", "newick")

    #  load gene df
    gdf = pd.read_csv(f"res/{hs}/gbk_annotations.csv")

    # load graph
    pan = pp.Pangraph.load_json(f"res/{hs}/joint_graph.json")

    # load tool annotations
    tdf = pd.read_csv(f"res/{hs}/tool_annotations.csv")

    # load hochhauser hotspots
    hh_genes = pd.read_csv("data/hh_genes.csv")
    up = hh_genes["Gene symbol of upstream gene"]
    down = hh_genes["Gene symbol of downstream gene"]
    # exclude nan
    up = set(up[~pd.isna(up)])
    down = set(down[~pd.isna(down)])

    return tree, gdf, pan, tdf, up, down


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


def draw_gene(ax, b, e, y, strand, w=1, c="black", arr_l=100):
    if strand:
        a = max(b, e - arr_l)
        X = [b, a, e, a, b, b]
        Y = [w / 2, w / 2, 0, -w / 2, -w / 2, w / 2]
    else:
        a = min(b + arr_l, e)
        X = [e, a, b, a, e, e]
        Y = [-w / 2, -w / 2, 0, w / 2, w / 2, -w / 2]
    Y = np.array(Y) + y
    ax.plot(X, Y, color=c, lw=1)


def draw_genes(ax, gdf, hh_up, hh_down, strain_y):

    colors = {
        "CDS": "black",
        "tRNA": "orange",
        "ncRNA": "green",
    }

    for i, row in gdf.iterrows():
        b, e, strand = row.start, row.end, row.strand
        kind, iso, gene = row.kind, row.iso, row.gene
        if kind == "gene":
            continue
        y = strain_y[iso]
        if gene in hh_up:
            c = "blue"
        elif gene in hh_down:
            c = "red"
        else:
            c = colors[kind]
        draw_gene(ax, b, e, y, strand, c=c, w=0.7)
        # add gene name if not NaN
        if not pd.isna(gene):
            gx = (b + e) / 2
            ax.text(gx, y, gene, fontsize=3, va="center", ha="left")


def mark_tool_ann(ax, tdf, strain_y):
    colors = {
        "ISEScan": "C0",
        "defensefinder": "C2",
        "integronfinder": "C3",
        "genomad": "C4",
    }

    # draw annotations as barh in the background
    for i, row in tdf.iterrows():
        y = strain_y[row.iso]

        ax.barh(
            y,
            row.end - row.start,
            left=row.start,
            color=colors[row.kind],
            alpha=0.5,
            height=0.5,
            zorder=-1,
        )


def mark_blocks(ax, strain_y, pan, block_colors):

    ax = axs[1]
    for iso, y in strain_y.items():
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
            if not strand:
                edgecolor = "k"
            ax.barh(
                y,
                x1 - x0,
                left=x0,
                height=0.8,
                color=color,
                edgecolor=edgecolor,
                linewidth=0.5,
                zorder=-1,
            )


if __name__ == "__main__":
    args = parse_args()
    hs = args.hs
    # hs = "CIRMBUYJFK_f__CWCCKOQCWZ_r"

    fig_fld = pathlib.Path(f"res/{hs}")
    fig_fld.mkdir(parents=True, exist_ok=True)

    tree, gdf, pan, tdf, hh_up, hh_down = load_data(hs)
    strain_y = {l.name: i + 1 for i, l in enumerate(tree.get_terminals())}

    gen = color_generator()
    block_colors = defaultdict(lambda: next(gen))

    # assign colors to first/last blocks
    strains = pan.strains()
    Bs = pan.paths[strains[0]].block_ids
    block_colors[Bs[0]] = left_clr
    block_colors[Bs[-1]] = right_clr

    fig, axs = plt.subplots(
        1, 2, figsize=(24, 20), sharey=True, gridspec_kw={"width_ratios": [1, 10]}
    )

    ax = axs[0]
    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: ""
    )
    despine(ax)
    ax.set_xlabel("branch length")
    ax.set_ylabel("isolates")

    ax = axs[1]
    draw_genes(ax, gdf, hh_up, hh_down, strain_y)
    mark_tool_ann(ax, tdf, strain_y)
    ax.set_xlabel("junction position (bp)")
    ax.set_xlim(0, max([pan.paths[iso].block_positions[-1] for iso in pan.strains()]))

    despine(ax)
    plt.tight_layout()
    plt.savefig(fig_fld / "gene_1.pdf")
    plt.savefig(fig_fld / "gene_1.png", dpi=200)
    show()

    fig, axs = plt.subplots(
        1, 2, figsize=(24, 20), sharey=True, gridspec_kw={"width_ratios": [1, 10]}
    )

    ax = axs[0]
    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: ""
    )
    despine(ax)

    ax = axs[1]
    draw_genes(ax, gdf, hh_up, hh_down, strain_y)
    mark_blocks(ax, strain_y, pan, block_colors)
    ax.set_xlabel("junction position (bp)")
    ax.set_xlim(0, max([pan.paths[iso].block_positions[-1] for iso in pan.strains()]))

    despine(ax)
    plt.tight_layout()
    plt.savefig(fig_fld / "gene_2.pdf")
    plt.savefig(fig_fld / "gene_2.png", dpi=200)
    show()
