import argparse
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser(description="Produce figures on distance matrices")
    parser.add_argument("--tree", type=str)
    parser.add_argument("--dist_df", type=str)
    parser.add_argument("--fig_fld", type=str)
    return parser.parse_args()


def to_dist_mat(df, leaf_order, value):
    D = df.pivot(index="si", columns="sj", values=value)
    D = D[leaf_order]
    D = D.loc[leaf_order]
    return D


def plot_dist_mat(D, tree, title, savefig):
    N = D.shape[0]
    fig, axs = plt.subplots(
        1,
        2,
        figsize=(3 + N * 0.2, 1 + N * 0.2),
        gridspec_kw={"width_ratios": [1, 4]},
    )

    # plot tree on left axis
    ax = axs[0]
    Phylo.draw(
        tree,
        axes=ax,
        show_confidence=False,
        label_func=lambda x: "",
        do_show=False,
    )
    ax.set_axis_off()

    # plot distance matrix on right axis
    ax = axs[1]
    ax.set_title(title)
    sns.heatmap(D, cmap="viridis", square=False, ax=ax, cbar_kws={"shrink": 0.3})
    plt.tight_layout()
    plt.savefig(savefig)
    plt.close()


values_to_titles = {
    # "si",
    # "sj",
    "core_div_naive": ("unfiltered core-genome divergence", "core_div_naive"),
    "mash_dist": ("mash distance", "mash_dist"),
    "private seq. (bp)": ("private seq. (bp)", "private_seq"),
    "shared seq. (bp)": ("shared seq. (bp)", "shared_seq"),
    # "n. breakpoints": ("n. breakpoints", "n_breakpoints"),
    # "part. entropy": ("partition entropy", "part_entropy"),
    "n. blocks": ("n. pairwise blocks", "n_pairwise_blocks"),
    "core_div_filtered": ("filtered core-genome divergence", "core_div_filtered"),
    "edge_PA": ("edge P/A distance", "edge_PA"),
    # "edge_PA_reduced": ("projection edge P/A distance", "edge_PA_reduced"),
    "edge_sharing": ("n. shared edges", "edge_sharing"),
    "block_PA": ("block P/A distance", "block_PA"),
    "block_sharing": ("n. shared blocks", "block_sharing"),
    # "acc_block_PA": ("accessory blocks P/A distance", "acc_block_PA"),
}

distributions = [
    ("private seq. (bp)", "mash_dist"),
    ("private seq. (bp)", "core_div_filtered"),
    ("private seq. (bp)", "n. blocks"),
    ("private seq. (bp)", "block_PA"),
    ("core_div_filtered", "core_div_naive"),
    ("core_div_filtered", "n. blocks"),
    ("core_div_filtered", "block_PA"),
    ("core_div_filtered", "edge_PA"),
    ("core_div_filtered", "private seq. (bp)"),
    ("core_div_filtered", "shared seq. (bp)"),
]

if __name__ == "__main__":
    args = parse_args()

    # save folder
    svfld = pathlib.Path(args.fig_fld)
    svfld.mkdir(exist_ok=True, parents=True)

    # load tree
    tree = Phylo.read(args.tree, format="newick")
    tree.root_at_midpoint()
    tree.ladderize()
    leaf_order = [l.name for l in tree.get_terminals()]

    # load distance matrix
    df = pd.read_csv(args.dist_df)

    # plot distance matrices
    for value, (title, sv_title) in values_to_titles.items():
        D = to_dist_mat(df, leaf_order, value)
        svname = svfld / f"{sv_title}.png"
        plot_dist_mat(D, tree, title, savefig=svname)

    # distributions
    L = len(distributions)
    nY = (L // 2) + int(L % 2)
    fig, axs = plt.subplots(nY, 2, figsize=(10, nY * 4))
    for n, (xv, yv) in enumerate(distributions):
        ax = axs[n // 2, n % 2]
        sns.histplot(data=df, x=xv, y=yv, bins=50, ax=ax)
        xlab, xtitle = values_to_titles[xv]
        ylab, ytitle = values_to_titles[yv]
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    plt.tight_layout()
    plt.savefig(svfld / "distributions.png")
    plt.close()
