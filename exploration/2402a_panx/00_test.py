# %%
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
from Bio import SeqIO, AlignIO, Phylo
import pathlib

fig_fld = pathlib.Path("figs/f00")
fig_fld.mkdir(exist_ok=True, parents=True)

# %%

# load alignment
aln = AlignIO.read(
    "../../data/panX/data/ST131_ABC/geneCluster/genePresence.aln", "fasta"
)

# load tree:
pangraph_tree_file = (
    "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
)
tree = Phylo.read(pangraph_tree_file, "newick")
strains = [l.name for l in tree.get_terminals()]

# load gene clusters info
gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"

with open(gene_cluster_json, "r") as f:
    gdf = pd.DataFrame(json.load(f))
gdf["divers"] = gdf["divers"].astype(float)
gdf["count"] = gdf["count"].astype(int)
gdf["geneLen"] = gdf["geneLen"].astype(int)
gdf["event"] = gdf["event"].astype(int)


# %%
M = {rec.id: list(rec) for rec in aln}
M = np.array([M[s] for s in strains])
N, L = M.shape

# convert to number
M = (M == "1").astype(int)
M
# %%
fig, axs = plt.subplots(
    2,
    2,
    figsize=(10, 10),
    gridspec_kw={
        "width_ratios": [0.15, 1],
        "height_ratios": [1, 10],
        # "hspace": 0.1,
        "wspace": 0.0,
    },
    # sharex="col",
    # sharey="row",
)

ax = axs[1, 0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
ax.set_ylabel("isolates")
sns.despine(ax=ax)

ax = axs[1, 1]
ax.matshow(M, cmap="Greys", aspect="auto")
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel("gene presence/absence")


ax = axs[0, 1]
ax.plot(M.sum(axis=0) / N, "k-")
ax.set_ylabel("gene freq.")
ax.set_xlim(0, L)
sns.despine(ax=ax)

ax = axs[0, 0]
ax.axis("off")

plt.tight_layout()
plt.savefig(fig_fld / "PA_matrix.png")
plt.show()

# %%
# compare trees:
panx_tree_file = "../../data/panX/data/ST131_ABC/vis/strain_tree.nwk"
pg_unfiltered_tree_file = "../../results/ST131_ABC/pangraph/asm20-100-5-coretree.nwk"

import subprocess

cmd = f"treeknit {panx_tree_file} {pangraph_tree_file} -o data/f00/tk_test --better-MCCs --auspice-view"
subprocess.run(cmd, shell=True)

cmd = f"treeknit {panx_tree_file} {pg_unfiltered_tree_file} -o data/f00/tk_test_unf --better-MCCs --auspice-view"
subprocess.run(cmd, shell=True)

# %%

fig, axs = plt.subplots(1, 2, figsize=(9, 4))

ax = axs[0]
sns.histplot(gdf, x="count", ax=ax, hue="dupli", bins=20, element="step")
ax.set_xlabel("n. isolates")
ax.set_ylabel("n. gene clusters")
sns.move_legend(ax, "upper center", title="duplicated")

ax = axs[1]
sns.histplot(
    gdf,
    x="divers",
    ax=ax,
    bins=[0] + list(np.logspace(-4, 0, 50)),
    element="step",
    cumulative=True,
)
ax.set_xscale("symlog", linthresh=0.001)
ax.set_xlim(0, 1)
ax.set_xlabel("diversity")
ax.set_ylabel("n. gene clusters")

plt.tight_layout()
plt.savefig(fig_fld / "gene_cluster_stats.png")
plt.show()

# %%

pan_file = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pan_file)
bdf = pan.to_blockstats_df()
# %%

fig, ax = plt.subplots(1, 1, figsize=(10, 4))

N = len(strains)
bins = np.arange(0, N + 1, 1) + 0.5

ax.hist(
    bdf["n. strains"],
    weights=bdf["len"],
    bins=bins,
    color="k",
    alpha=0.5,
    cumulative=True,
    density=True,
    histtype="step",
    label="pangraph",
)

ax.hist(
    gdf["count"],
    weights=gdf["geneLen"],
    bins=bins,
    color="r",
    alpha=0.5,
    cumulative=True,
    density=True,
    histtype="step",
    label="panx",
)
ax.set_xlim(0, N + 1)
ax.set_xlabel("n. strains")
ax.set_ylabel("pangenome fraction")

# add N to xticks
xticks = list(ax.get_xticks())
xticks.append(N)
ax.set_xticks(xticks)
ax.set_xlim(0, N + 1)

ax.legend(loc="upper left")
sns.despine()


plt.tight_layout()
plt.savefig(fig_fld / "panx_pangraph_pangenome_freq.png")
plt.savefig(fig_fld / "panx_pangraph_pangenome_freq.svg")
plt.savefig(fig_fld / "panx_pangraph_pangenome_freq.pdf")
plt.show()
# %%

pg_l = bdf["len"].sum()
print(f"pangraph pangenome size = {pg_l:_}")

panx_l = gdf["geneLen"].sum()
print(f"panx pangenome size = {panx_l:_}")

print("total number of gene clusters:")
print(len(gdf))

print("number of gene clusters:")
print(gdf["dupli"].value_counts())

print("number of singletons:")
mask = gdf["count"] == 1
mask &= gdf["dupli"] == "no"
print(mask.sum())

# %%
bins = [0, 1, 222 * 0.05, 222 * 0.30, 222 * 0.70, 222 * 0.95, 221, 222]

# split in count categories
mask = gdf["dupli"] == "no"
ct = pd.cut(gdf["count"][mask], bins=bins).value_counts()
pd.DataFrame([ct, ct / ct.sum() * 100]).T


# %%
import matplotlib as mpl

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
sns.histplot(
    gdf[gdf["dupli"] == "no"],
    y="event",
    x="count",
    discrete=True,
    ax=ax,
    cbar=True,
    cbar_kws={"label": "n. gene clusters", "extend": "max"},
    vmin=0,
    vmax=50,
)
plt.plot([0, N / 2, N], [0, N / 2, 0], ":", c="gray", alpha=0.3)
ax.set_xlabel("n. isolates")
ax.set_ylabel("n. events")
sns.despine()
ax.set_xlim(0, N)
ax.set_ylim(0, N / 2)
plt.tight_layout()
plt.savefig(fig_fld / "gene_cluster_event_count.png")
plt.show()
# %%
