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

distributions = {
    "private seq. (bp)": (
        "mash_dist",
        "core_div_filtered",
        "n. blocks",
        "block_PA",
    ),
    "core_div_filtered": (
        "core_div_naive",
        "n. blocks",
        "block_PA",
        "edge_PA",
        "private seq. (bp)",
        "shared seq. (bp)",
    ),
}


def plot_distribution(df, x_value, Y_values, values_to_titles, svname):
    L = len(Y_values)
    nY = (L // 2) + int(L % 2)
    fig, axs = plt.subplots(nY, 2, figsize=(8, nY * 4))
    mask = df["si"] > df["sj"]
    sdf = df[mask]
    for n, yv in enumerate(Y_values):
        ax = axs[n // 2, n % 2]
        sns.histplot(data=sdf, x=x_value, y=yv, bins=50, ax=ax)
        xlab = values_to_titles[x_value][0]
        ylab = values_to_titles[yv][0]
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
    plt.tight_layout()
    plt.savefig(svname)
    plt.close()


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
    for x_value, Y_values in distributions.items():
        x_title = values_to_titles[x_value][1]
        svname = svfld / f"distributions_{x_title}.png"
        plot_distribution(df, x_value, Y_values, values_to_titles, svname)
