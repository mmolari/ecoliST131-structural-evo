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
# E_dist = dist_df.reset_index().pivot(index="si", columns="sj", values="edge_sharing")
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

# cmap = plt.get_cmap("Greys_r")
cmap = plt.get_cmap("viridis_r")

ax = axs[0]
g = ax.imshow(T_dist, vmin=0, cmap=cmap)
ax.set_title("Core genome divergence")
plt.colorbar(g, ax=ax, shrink=0.6)

ax = axs[1]
g = ax.imshow(B_dist, vmin=0, cmap=cmap)
ax.set_title("Block presence/absence distance")
plt.colorbar(g, ax=ax, shrink=0.6)

ax = axs[2]
g = ax.imshow(E_dist, vmin=0, cmap=cmap)
ax.set_title("Edge presence/absence distance")
plt.colorbar(g, ax=ax, shrink=0.6)

for ax in axs:
    ax.set_xlabel("Strain")
    ax.set_ylabel("Strain")

plt.tight_layout()
svfig("distance_matrices")
plt.show()

# %%

fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharex=True)

ax = axs[0]
# ax.hist2d(T_dist.values.flatten(), E_dist.values.flatten(), bins=100, cmap="Blues")
ax.scatter(T_dist.values.flatten(), E_dist.values.flatten(), s=1, alpha=0.3)
ax.set_xlabel("Core genome divergence")
ax.set_ylabel("Edge presence/absence distance")

ax = axs[1]
# ax.hist2d(T_dist.values.flatten(), B_dist.values.flatten(), bins=100, cmap="Blues")
ax.scatter(T_dist.values.flatten(), B_dist.values.flatten(), s=1, alpha=0.3)
ax.set_xlabel("Core genome divergence")
ax.set_ylabel("Block presence/absence distance")

plt.tight_layout()
# svfig("scatterplot_snps_vs_edges_and_blocks")
plt.show()

# %%

# triangular histogram plot with seaborn

import seaborn as sns

sns.set_theme(style="darkgrid")

df = dist_df[
    [
        "core_div_filtered",
        "block_PA",
        "edge_PA",
        "edge_PA_reduced",
        "block_sharing",
        "edge_sharing",
    ]
]
g = sns.PairGrid(df, diag_sharey=False)
g.map_lower(sns.histplot)
g.map_diag(sns.histplot)
for i, j in zip(*np.triu_indices_from(g.axes, k=1)):
    g.axes[i, j].set_visible(False)
plt.show()


df = dist_df[
    [
        "core_div_filtered",
        "block_PA",
        "edge_PA",
        # "edge_PA_reduced",
        # "block_sharing",
        # "edge_sharing",
    ]
]
g = sns.PairGrid(df, diag_sharey=False)
g.map_lower(sns.histplot)
g.map_diag(sns.histplot)
for i, j in zip(*np.triu_indices_from(g.axes, k=1)):
    g.axes[i, j].set_visible(False)
svfig("pairgrid_distances")
plt.show()


# %%

import tanglegram as tg

# Plot tanglegram
fig = tg.plot(
    B_dist,
    E_dist,
    # sort="permutations",
    sort="step2side",
    figsize=(12, 15),
    # dend_kwargs={"color_threshold": 200},
)
plt.show()

# %%

# Plot tanglegram
fig = tg.plot(
    T_dist,
    B_dist,
    # sort="permutations",
    sort="step2side",
    figsize=(12, 15),
    # dend_kwargs={"color_threshold": 200},
)
plt.show()

# %%


# Plot tanglegram
fig = tg.plot(
    T_dist,
    E_dist,
    # sort="permutations",
    sort="step2side",
    figsize=(12, 15),
    # dend_kwargs={"color_threshold": 200},
)
plt.show()

# %%
