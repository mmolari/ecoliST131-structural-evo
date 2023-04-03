# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import numpy as np
from Bio import Phylo
import itertools as itt
import json

# %%

fig_fld = pathlib.Path("figs")

tree_old = "../../results/ST131/pangraph/asm20-100-5-coretree.nwk"
tree_new = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
dist_file = "../../results/ST131/distances/summary-asm20-100-5.csv"

info_old = "../../results/ST131/pangraph/asm20-100-5-alignment/corealignment_info.json"
info_new = "../../results/ST131/pangraph/asm20-100-5-alignment/filtered_corealignment_info.json"
# %%

with open(info_old, "r") as f:
    info = json.load(f)
factor_old = info["n. consensus"] + info["n. snps"]

with open(info_new, "r") as f:
    info = json.load(f)
factor_new = info["polished aln consensus"] + info["polished aln snps"]
# %%


def tree_to_distmat(tree_file, label, factor):
    tree = Phylo.read(tree_file, format="newick")
    leaves = tree.get_terminals()
    df = []
    for a, b in itt.product(leaves, leaves):
        df.append({"si": a.name, "sj": b.name, label: tree.distance(a, b) * factor})
    return pd.DataFrame(df).set_index(["si", "sj"])


dfo = tree_to_distmat(tree_old, "treedist_naive", factor=factor_old)
dfn = tree_to_distmat(tree_new, "treedist_refined", factor=factor_new)

# %%

ddf = pd.read_csv(dist_file).set_index(["si", "sj"])

df = pd.concat([dfo, dfn, ddf], axis=1)
df["core_div_naive"] *= factor_old
df["core_div_filtered"] *= factor_new

mask = df.index.get_level_values(0) > df.index.get_level_values(1)
df = df[mask]
# %%


def compare_plot(df, savename=None, max_div=None):
    fig, axs = plt.subplots(3, 2, figsize=(9, 12))

    if max_div is None:
        sdf = df
    else:
        mask = df["core_div_naive"] <= max_div
        sdf = df[mask]

    ax = axs[0, 0]
    M = sdf["core_div_naive"].max()
    ax.plot([0, M], [0, M], "k--", zorder=-1)
    ax.scatter(sdf["core_div_naive"], sdf["treedist_naive"], alpha=0.1)
    ax.set_xlabel("corealn snps (naive)")
    ax.set_ylabel("tree distance (naive)")

    ax = axs[0, 1]
    M = sdf["treedist_refined"].max()
    ax.plot([0, M], [0, M], "k--", zorder=-1)
    ax.scatter(sdf["core_div_filtered"], sdf["treedist_refined"], alpha=0.1)
    ax.set_ylabel("tree distance (polished)")

    ax = axs[1, 0]
    ax.hist2d(sdf["core_div_naive"], sdf["treedist_naive"], bins=50)
    ax.set_ylabel("tree distance (naive)")

    ax = axs[1, 1]
    ax.hist2d(sdf["core_div_filtered"], sdf["treedist_refined"], bins=50)
    ax.set_ylabel("tree distance (polished)")

    ax = axs[2, 0]
    ax.hist2d(sdf["core_div_naive"], sdf["mash_dist"], bins=50)
    ax.set_ylabel("mash distance")

    ax = axs[2, 1]
    ax.hist2d(sdf["core_div_filtered"], sdf["mash_dist"], bins=50)
    ax.set_ylabel("mash distance")

    for ax in axs[:, 0]:
        ax.set_xlabel("corealn snps (naive)")
    for ax in axs[:, 1]:
        ax.set_xlabel("corealn snps (polished)")

    sns.despine(fig)
    plt.tight_layout()
    if savename is not None:
        plt.savefig(fig_fld / savename)
    plt.show()


compare_plot(df, savename="benchmark_refinement.pdf")
compare_plot(df, max_div=500, savename="benchmark_refinement_masked.pdf")
# %%

sns.histplot(data=df, x="core_div_naive", y="core_div_filtered")
sns.despine()
plt.xlabel("core snps (naive)")
plt.ylabel("core snps (filtered)")
plt.tight_layout()
plt.savefig(fig_fld / "naive_vs_refined_dist.pdf")
plt.show()

# %%

from Bio import AlignIO

# %%
aln = AlignIO.read(
    "../../results/ST131/pangraph/asm20-100-5-alignment/filtered_corealignment.fa",
    format="fasta",
)
aln = np.array(aln)
# %%
N, L = aln.shape
# %%
Ms = []
for i, j in itt.product(range(N), range(N)):
    if i >= j:
        continue
    m = np.sum(aln[i, :] != aln[j, :])
    Ms.append(m)

# %%
plt.hist(Ms, bins=100)
plt.hist(df["core_div_filtered"], bins=100, histtype="step")
# %%
