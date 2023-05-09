# %%

import pathlib

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from Bio import Phylo


fig_p = pathlib.Path("fig")
fig_p.mkdir(exist_ok=True)


def svfig(svname):
    for k in ["pdf", "png", "svg"]:
        plt.savefig(fig_p / f"{svname}.{k}", dpi=300)


# %%

tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, format="newick")
tree.ladderize()

# %%
dist_file = "../../results/ST131/distances/summary-asm20-100-5.csv"
dist_df = pd.read_csv(dist_file, index_col=[0, 1])
T_dist = dist_df.reset_index().pivot(
    index="si", columns="sj", values="core_div_filtered"
)
E_dist = dist_df.reset_index().pivot(index="si", columns="sj", values="edge_PA")
# E_dist = dist_df.reset_index().pivot(index="si", columns="sj", values="edge_PA_reduced")
B_dist = dist_df.reset_index().pivot(index="si", columns="sj", values="block_PA")

# %%
leaves_order = [l.name for l in tree.get_terminals()]
B_dist = B_dist.reindex(leaves_order, axis=0).reindex(leaves_order, axis=1)
E_dist = E_dist.reindex(leaves_order, axis=0).reindex(leaves_order, axis=1)
T_dist = T_dist.reindex(leaves_order, axis=0).reindex(leaves_order, axis=1)

# %%
# plot
fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)

ax = axs[0]
g = ax.imshow(B_dist, cmap="Greys", vmin=0)
ax.set_title("Block presence/absence distance")
plt.colorbar(g, ax=ax, shrink=0.6)

ax = axs[1]
g = ax.imshow(E_dist, cmap="Greys", vmin=0)
ax.set_title("Edge presence/absence distance")
plt.colorbar(g, ax=ax, shrink=0.6)

ax = axs[2]
g = ax.imshow(T_dist, cmap="Greys", vmin=0)
ax.set_title("Core genome divergence")
plt.colorbar(g, ax=ax, shrink=0.6)

for ax in axs:
    ax.set_xlabel("Strain")
    ax.set_ylabel("Strain")

plt.tight_layout()
svfig("distance_matrices")
plt.show()

# %%
