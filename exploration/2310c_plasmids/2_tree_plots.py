# %%
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from Bio import Phylo
import argparse

df_chr_file = "../../results/ST131_full/resistance/ncbi_summary.txt"
df_pls_file = "../../results/ST131_full/plasmids/resistance/ncbi_summary_chromosome.csv"
tree_file = "../../results/ST131_full/pangraph/asm20-100-5-filtered-coretree.nwk"

# %%


def parse_args():
    parser = argparse.ArgumentParser(description="Plot metadata")
    parser.add_argument("--coregenome_tree", type=str)
    parser.add_argument("--chr_df", type=str)
    parser.add_argument("--pls_df", type=str)
    parser.add_argument("--id_threshold", type=float)
    parser.add_argument("--fig", type=str)
    return parser.parse_args()


def parse_element(x, thr):
    """
    Assign `.` = 0
    and `100.00;100;00;100.00` = 3
    """
    if x == ".":
        return 0
    x = str(x)
    els = x.split(";")
    ct = 0
    for el in els:
        if float(el) > thr * 100.0:
            ct += 1
    return ct


def load_res_df(fname, thr):
    df = pd.read_csv(fname, sep="\t", index_col=0)
    # transform filename to accession number
    df.index = df.index.map(lambda x: pathlib.Path(x).stem)
    df.drop(columns="NUM_FOUND", inplace=True)
    # parse elements
    df = df.applymap(lambda x: parse_element(x, thr))
    return df


def plot_resistance_vs_tree(tree, df_chr, df_pls):
    I = list([l.name for l in tree.get_terminals()])
    R = np.unique(list(df_chr.columns) + list(df_pls.columns))
    X, Y = len(R), len(I)

    iso_y = {l: n + 1 for n, l in enumerate(I)}

    fig, axs = plt.subplots(
        1,
        2,
        figsize=(0.2 * X + 3, 1 + 0.08 * Y),
        sharey=True,
        gridspec_kw={"width_ratios": [1, X * 0.10]},
    )

    # plot tree
    ax = axs[0]
    Phylo.draw(tree, lambda x: "", do_show=False, axes=ax)

    # plot P/A matrix
    ax = axs[1]
    xticks = {}
    x = 0
    for r in R:
        xticks[x] = r
        for iso in I:
            y = iso_y[iso]
            # chromosome
            if r in df_chr.columns and iso in df_chr.index:
                ct = df_chr.loc[iso, r]
            else:
                ct = 0
            if ct > 0:
                ax.scatter(x, y, color="blue", marker="3")

            # plasmid
            if r in df_pls.columns and iso in df_pls.index:
                ct = df_pls.loc[iso, r]
            else:
                ct = 0
            if ct > 0:
                ax.scatter(x, y, color="red", marker="4")
        ax.axvline(x, c="lightgray", zorder=-1, alpha=0.3)
        x += 1
    ax.set_xticks(list(xticks.keys()))
    ax.set_xticklabels(list(xticks.values()), rotation=90)
    ax.set_title("n. resistance genes")

    for ax in axs:
        for n in range(Y):
            ax.axhline(n + 1, c="lightgray", zorder=-1, alpha=0.3)
        ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    return fig, axs


# %%

df_chr = load_res_df(df_chr_file, 0.95)
df_pls = pd.read_csv(df_pls_file, index_col=0)
tree = Phylo.read(tree_file, "newick")

# %%


plot_resistance_vs_tree(tree, df_chr, df_pls)

# %%
