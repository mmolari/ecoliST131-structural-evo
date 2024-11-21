# %%

import json
import treetime
import pathlib

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib as mpl

from Bio import Phylo, AlignIO
from collections import Counter

res_pth = pathlib.Path("../../results/ST131_ABC/pangraph")


def parse_args():
    return {
        "tree": res_pth / "asm20-100-5-coretree.nwk",
        "aln": res_pth / "asm20-100-5-alignment/corealignment.fa",
        "aln_info": res_pth / "asm20-100-5-alignment/corealignment_info.json",
        "filt_tree": res_pth / "asm20-100-5-filtered-coretree.nwk",
        "filt_aln": res_pth / "asm20-100-5-alignment/filtered_corealignment.fa",
        "filt_aln_info": res_pth
        / "asm20-100-5-alignment/filtered_corealignment_info_size.json",
    }


def treetime_infer(res_dict):
    tree = Phylo.read(res_dict["tree_file"], format="newick")
    aln_file = res_dict["aln_f"]
    aln_L = res_dict["aln_L"]

    # instantiate treetime
    myTree = treetime.TreeAnc(
        gtr="Jukes-Cantor", tree=tree, aln=aln_file, verbose=0, seq_len=aln_L
    )
    myTree.tree.root.branch_length = 0.0
    myTree.infer_ancestral_sequences(prune_short=True)

    return myTree


def count_mutations_df(tt_tree):
    mut_count = Counter()
    for node in tt_tree.tree.find_clades():
        muts = node.mutations
        mut_count.update([m[1] for m in muts])

    df = pd.Series(mut_count, name="n. muts").reset_index()
    df = df.sort_values("index").reset_index()
    return df


def aln_biallelic_idxs(aln_file):
    def __is_biallelic(row):
        vals, cts = np.unique(row, return_counts=True)
        return (len(vals) == 2) and (cts.min() > 1)

    aln = AlignIO.read(aln_file, format="fasta")
    aln = np.array(aln)
    mask = np.apply_along_axis(__is_biallelic, axis=0, arr=aln)
    idxs = np.argwhere(mask).flatten()
    return idxs


def aln_terminal_idxs(aln_file):
    def __is_terminal(row):
        vals, cts = np.unique(row, return_counts=True)
        return (len(vals) == 2) and (cts.min() == 1)

    aln = AlignIO.read(aln_file, format="fasta")
    aln = np.array(aln)
    mask = np.apply_along_axis(__is_terminal, axis=0, arr=aln)
    idxs = np.argwhere(mask).flatten()
    return idxs


def mut_df(res_b, res_a):
    res_b["TT"] = treetime_infer(res_b)
    res_a["TT"] = treetime_infer(res_a)

    dfb = count_mutations_df(res_b["TT"])
    dfa = count_mutations_df(res_a["TT"])

    bidxsa = aln_biallelic_idxs(res_a["aln_f"])
    bidxsb = aln_biallelic_idxs(res_b["aln_f"])

    tidxsa = aln_terminal_idxs(res_a["aln_f"])
    tidxsb = aln_terminal_idxs(res_b["aln_f"])

    dfb["aln"] = "raw"
    dfa["aln"] = "filtered"

    dfa["phylogenetic info"] = "multi-allelic"
    dfb["phylogenetic info"] = "multi-allelic"

    mask_a = dfa["index"].isin(bidxsa).to_numpy()
    dfa.loc[mask_a, "phylogenetic info"] = "bi-allelic"
    mask_b = dfb["index"].isin(bidxsb).to_numpy()
    dfb.loc[mask_b, "phylogenetic info"] = "bi-allelic"

    mask_a = dfa["index"].isin(tidxsa).to_numpy()
    dfa.loc[mask_a, "phylogenetic info"] = "terminal branch"
    mask_b = dfb["index"].isin(tidxsb).to_numpy()
    dfb.loc[mask_b, "phylogenetic info"] = "terminal branch"

    df = pd.concat([dfa, dfb], ignore_index=True).drop(columns="level_0")
    df = df.rename(columns={"index": "position"})
    df["aln"] = pd.Categorical(df["aln"], categories=["raw", "filtered"])
    df["homoplasy"] = df["n. muts"] > 1
    return df


# def fig_hist(df, svname):
# g = sns.FacetGrid(
#     datas,
#     col="aln",
#     row="phylogenetic info",
#     sharex="col",
#     sharey="row",
#     height=3,
#     aspect=1.5,
# )

# g.map_dataframe(
#     sns.histplot,
#     x="position",
#     element="step",
#     binwidth=100,
#     hue="homoplasy",
#     hue_order=[False, True],
# )

# def add_multiallelic_frac(data, color, *args, **kwargs):
#     ax = plt.gca()
#     t = ax.transAxes
#     frac = data["homoplasy"].mean()
#     ns = (~data["homoplasy"]).sum()
#     nm = data["homoplasy"].sum()
#     ax.text(
#         0.05,
#         0.9,
#         f"homopl. frac. = {frac*100:.1f} %",
#         color="darkorchid",
#         transform=t,
#     )
#     ax.text(0.55, 0.9, f"single mut. (n = {ns})", color="C0", transform=t)
#     ax.text(0.55, 0.83, f" homoplasy (n = {nm})", color="C1", transform=t)
#     ax.set_ylim(top=ax.get_ylim()[1] * 1.1)

# g.map_dataframe(add_multiallelic_frac)
# plt.tight_layout()
# plt.savefig(svname)
# plt.close()

# %%

args = parse_args()

# load results before filtering
with open(args["aln_info"], "r") as f:
    info = json.load(f)
aln_L = info["n. consensus"] + info["n. snps"]

res_b = {
    "tree_file": str(args["tree"]),
    "aln_f": str(args["aln"]),
    "aln_L": aln_L,
}

# load results after filtering
with open(args["filt_aln_info"], "r") as f:
    info = json.load(f)
aln_L = info["polished aln size"]

res_a = {
    "tree_file": str(args["filt_tree"]),
    "aln_f": str(args["filt_aln"]),
    "aln_L": aln_L,
}

# infer ancestral sequences, count mutations and find homoplasies
df = mut_df(res_b, res_a)
# %%


def fig_hist(df):
    mask = df["phylogenetic info"] == "bi-allelic"
    sdf = df[mask]
    sdf

    fig, axs = plt.subplots(1, 2, figsize=(9, 4), sharey=True)

    for ax, tp, leg, bw, wh in [
        (axs[0], "raw", False, 300, "pre"),
        (axs[1], "filtered", True, 100, "post"),
    ]:
        mask = sdf["aln"] == tp
        sns.histplot(
            sdf[mask],
            x="position",
            element="step",
            binwidth=bw,
            weights=1 / bw,
            hue="homoplasy",
            multiple="stack",
            hue_order=[False, True],
            legend=leg,
            ax=ax,
        )
        ax.set_ylabel("fraction of bi-allelic non-singleton substitutions")
        ax.set_xlabel("restricted alignment position")
        ax.set_title(f"{wh} recombination filter")
        h_frac = sdf[mask]["homoplasy"].mean()
        ax.text(
            0.05,
            0.95,
            f"fraction of homoplasies = {h_frac*100:.1f} %",
            transform=ax.transAxes,
            color="sienna",
            fontsize=12,
        )

    ax.set_ylim(0, 1)

    sns.despine()

    return fig, axs


fig, axs = fig_hist(df)
plt.tight_layout()
plt.savefig("figs/f02/homoplasies.pdf")
plt.show()

# %%
print(df["aln"].value_counts().sort_index())
print(df[["aln", "phylogenetic info"]].value_counts().sort_index())
print(
    df[["aln", "phylogenetic info"]].value_counts().sort_index()
    / df["aln"].value_counts()
)


# %%
print(df[["aln", "phylogenetic info", "homoplasy"]].value_counts().sort_index())

# %%
