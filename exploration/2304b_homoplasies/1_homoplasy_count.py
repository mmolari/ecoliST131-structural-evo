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

fig_p = pathlib.Path("fig")
fig_p.mkdir(exist_ok=True)


def svfig(svname):
    for k in ["pdf", "png", "svg"]:
        plt.savefig(fig_p / f"{svname}.{k}", dpi=300)


res_f = pathlib.Path("../../results/ST131/pangraph")

# load original tree
tree_file = res_f / "asm20-100-5-coretree.nwk"
aln_file = res_f / "asm20-100-5-alignment/corealignment.fa"
aln_info_file = res_f / "asm20-100-5-alignment/corealignment_info.json"

with open(aln_info_file, "r") as f:
    info = json.load(f)
aln_L = info["n. consensus"] + info["n. snps"]

res_b = {
    "tree_file": str(tree_file),
    "aln_f": str(aln_file),
    "aln_L": aln_L,
}

# load filtered tree

tree_file = res_f / "asm20-100-5-filtered-coretree.nwk"
aln_file = res_f / "asm20-100-5-alignment/filtered_corealignment.fa"
aln_info_file = res_f / "asm20-100-5-alignment/filtered_corealignment_info.json"

with open(aln_info_file, "r") as f:
    info = json.load(f)
aln_L = info["polished aln size"]

res_a = {
    "tree_file": str(tree_file),
    "aln_f": str(aln_file),
    "aln_L": aln_L,
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


# %%

res_b["TT"] = treetime_infer(res_b)
res_a["TT"] = treetime_infer(res_a)

dfb = count_mutations_df(res_b["TT"])
dfa = count_mutations_df(res_a["TT"])

bidxsa = aln_biallelic_idxs(res_a["aln_f"])
bidxsb = aln_biallelic_idxs(res_b["aln_f"])

tidxsa = aln_terminal_idxs(res_a["aln_f"])
tidxsb = aln_terminal_idxs(res_b["aln_f"])
# %%
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
df["multimut"] = df["n. muts"] > 1

# %%

g = sns.FacetGrid(
    data=df,
    col="aln",
    row="phylogenetic info",
    # hue="multimut",
    sharex="col",
    sharey="row",
    height=3,
    aspect=1.5,
    # subplot_kws={"figsize": (8, 8)},
)
g.map_dataframe(
    sns.histplot,
    x="position",
    element="step",
    binwidth=100,
    hue="multimut",
    hue_order=[False, True],
)


def add_multiallelic_frac(data, color, *args, **kwargs):
    ax = plt.gca()
    t = ax.transAxes
    frac = data["multimut"].mean()
    ns = (~data["multimut"]).sum()
    nm = data["multimut"].sum()
    ax.text(
        0.05, 0.9, f"homopl. frac. = {frac*100:.1f} %", color="darkorchid", transform=t
    )
    ax.text(0.55, 0.9, f"single mut. (n = {ns})", color="C0", transform=t)
    ax.text(0.55, 0.83, f" homoplasy (n = {nm})", color="C1", transform=t)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.1)


g.map_dataframe(add_multiallelic_frac)
plt.tight_layout()
svfig("homoplasies")
plt.show()

# %%
# mask = (df["aln"] == "raw") & (df["phylogenetic info"] == "bi-allelic")
# mask &= df["multimut"]
# I = df[mask]["position"].to_numpy()
mask = (df["aln"] == "filtered") & (df["phylogenetic info"] == "bi-allelic")
mask &= df["multimut"]
I = df[mask]["position"].to_numpy()
# %%


def is_cons(row):
    vals, cts = np.unique(row, return_counts=True)
    cons_i = np.argmax(cts).flatten()
    return row == vals[cons_i]


def subaln(aln_f, idxs, str_order):
    aln = AlignIO.read(aln_f, format="fasta")
    strains = [a.name for a in aln]
    strain_idx = {s: n for n, s in enumerate(strains)}
    aln = np.array(aln)
    order = np.array([strain_idx[s] for s in str_order])
    saln = aln[order][:, idxs]
    return np.apply_along_axis(is_cons, axis=0, arr=saln)


# %%
tree = Phylo.read(res_a["tree_file"], format="newick")
tree.ladderize()
str_ord = [n.name for n in tree.get_terminals()]
saln = subaln(res_a["aln_f"], idxs=I, str_order=str_ord)

fig, axs = plt.subplots(1, 2, sharey=True, figsize=(5, 8))
ax = axs[0]
Phylo.draw(tree, label_func=lambda x: None, do_show=False, axes=ax)
ax.set_title("core genome tree")


ax = axs[1]
aw = np.argwhere(~saln)
ax.scatter(y=aw[:, 0] + 1, x=aw[:, 1], c=aw[:, 1], cmap="tab20")
ax.set_ylim(top=-1)
K = saln.shape[1]
ax.set_xticks(np.arange(K), minor=True)
ax.set_title("homoplasies")

for ax in axs:
    N = len(str_ord)
    ax.set_yticks(np.arange(N) + 1, minor=True)
    ax.grid(which="major", alpha=0.6)
    ax.grid(which="minor", alpha=0.3)

plt.tight_layout()
svfig("remaining_biallelic_homoplasies")
plt.show()
# %%

# color supported branches
mask = (df["aln"] == "filtered") & (df["phylogenetic info"] == "bi-allelic")
mask &= ~df["multimut"]
I = df[mask]["position"].to_numpy()

T = res_a["TT"].tree
Ls = [len([m for m in n.mutations if m[1] in I]) for n in T.get_nonterminals()]
cmap = plt.get_cmap("Blues")
norm = mpl.colors.LogNorm(vmin=1, vmax=max(Ls))
# cmap = plt.get_cmap("cool")
# norm = mpl.colors.BoundaryNorm(boundaries=[0,1,10,100,300], ncolors=256)

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
for n in T.get_nonterminals():
    k = len([m for m in n.mutations if m[1] in I])
    c = cmap(norm(k))
    n.color = mpl.colors.to_hex(c)
for n in T.get_terminals():
    n.color = "#888888"
T.root.color = "#888888"
Phylo.draw(T, label_func=lambda x: None, do_show=False, axes=ax)
mapp = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
plt.colorbar(
    mapp,
    label="n. supporting bi-allelic SNPs",
    shrink=0.8,
    extend="min",
)
sns.despine()
plt.tight_layout()
svfig("branch_support")
plt.show()
# %%
