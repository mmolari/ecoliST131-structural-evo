# %%
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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

cmd = f"treeknit {panx_tree_file} {pangraph_tree_file} -o tk_test --better-MCCs --auspice-view"
subprocess.run(cmd, shell=True)

cmd = f"treeknit {panx_tree_file} {pg_unfiltered_tree_file} -o tk_test_unf --better-MCCs --auspice-view"
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

# %%
