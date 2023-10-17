import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
import json
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Plot metadata")
    parser.add_argument("--coregenome_tree", type=str)
    parser.add_argument("--mlst_df", type=str)
    parser.add_argument("--plsm_json", type=str)
    parser.add_argument("--fig", type=str)
    return parser.parse_args()


def parse_tree(tree_file):
    tree = Phylo.read(tree_file, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    return tree


def parse_mlst_df(mlst_df_file):
    df = pd.read_csv(mlst_df_file, sep=",", index_col=0, dtype=str)
    df.fillna("-", inplace=True)
    df["typ"] = "F" + df["FII"] + ":A" + df["FIA"] + ":B" + df["FIB"]

    # remove untyped
    mask = df["typ"] == "F-:A-:B-"
    df = df[~mask]
    return df


def reduce_mlst_df(tree, mlst_df, n_top):
    I = [l.name for l in tree.get_terminals()]
    J = mlst_df["typ"].value_counts().index[:n_top]

    pa_df = []
    for iso in I:
        res = {"iso": iso}
        if not iso in pls_dict:
            pa_df.append(res)
            continue
        P = pls_dict[iso]
        for plsm in P:
            if not plsm in mlst_df.index:
                continue
            t = mlst_df.loc[plsm, "typ"]
            if t in J:
                res[t] = True
        pa_df.append(res)
    pa_df = pd.DataFrame(pa_df)
    pa_df.set_index("iso", inplace=True)
    pa_df.fillna(False, inplace=True)
    pa_df = pa_df[J]
    pa_df = pa_df.loc[I]
    return pa_df


def make_fig(tree, pa_df, fig_file):
    N = len(tree.get_terminals())
    fig, axs = plt.subplots(1, 2, figsize=(10, N * 0.05))
    ax = axs[0]
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
    ax.set_yticks([])
    # ax.set_xticks([])
    ax.set_title("ST131 core tree")

    ax = axs[1]
    sns.heatmap(pa_df, ax=ax, cbar=False, cmap="binary")
    ax.set_yticks([])
    ax.set_ylabel("Isolate")
    ax.set_xlabel("Plasmid type")
    ax.set_title("Plasmid types in ST131 isolates")
    plt.tight_layout()
    sns.despine()
    plt.savefig(fig_file)
    plt.close()


if __name__ == "__main__":
    args = parse_args()
    tree = parse_tree(args.coregenome_tree)
    mlst_df = parse_mlst_df(args.mlst_df)

    with open(args.plsm_json, "r") as f:
        pls_dict = json.load(f)

    pa_df = reduce_mlst_df(tree, mlst_df, 8)

    make_fig(tree, pa_df, args.fig)
