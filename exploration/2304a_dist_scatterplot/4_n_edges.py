# %%
import pathlib
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import Phylo

fig_fld = pathlib.Path("figs")
fig_fld.mkdir(exist_ok=True)


def svfig(svname):
    for suffix in ["png", "pdf"]:
        plt.savefig(
            fig_fld / f"{svname}.{suffix}",
            dpi=300,
            facecolor="white",
        )


tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, format="newick")
tree.ladderize()
str_ord = [x.name for x in tree.get_terminals()]

df_file = "../../results/ST131/distances/summary-asm20-100-5.csv"
adf = pd.read_csv(df_file)
mask = adf["si"] > adf["sj"]
df = adf[mask].copy()


df.rename(
    columns={
        "core_div_filtered": "core genome div.",
        "edge_PA_reduced": "edge_PA_reduced",
        "block_PA": "block P/A",
        "n. blocks": "n. pairwise blocks",
        "edge_PA": "edge P/A",
        "edge_PA_reduced": "edge P/A (reduced)",
    },
    inplace=True,
)


cols = [
    "si",
    "sj",
    "core_div_naive",
    "mash_dist",
    "private seq. (bp)",
    "shared seq. (bp)",
    "n. breakpoints",
    "part. entropy",
    "n. pairwise blocks",
    "core genome div.",
    "edge_PA",
    "edge_PA_reduced",
    "edge_sharing",
    "block P/A",
    "block_sharing",
]

# %%

variables = [
    "core genome div.",
    # "block P/A",
    "edge P/A",
    "edge P/A (reduced)",
    # "edge_sharing",
    # "n. pairwise blocks",
    # "block P/A",
    # "block_sharing",
    # "part. entropy",
    # "private seq. (bp)",
]

color_strain = "NZ_JAOSEJ010000001"
df["NZ_JAOSEJ010000001"] = (df["si"] == color_strain) | (df["sj"] == color_strain)


g = sns.PairGrid(
    df,
    # hue="NZ_JAOSEJ010000001",
    vars=variables,
    diag_sharey=False,
)
g.map_lower(sns.histplot)
g.map_diag(sns.histplot)
# g.map_upper(sns.kdeplot, fill=True)

# hide upper-diagonal axes
for i, j in zip(*np.triu_indices_from(g.axes, 1)):
    g.axes[i, j].set_visible(False)

plt.tight_layout()
svfig("6_edges_vs_core")
plt.show()


# %%


fig, axs = plt.subplots(1, 2, figsize=(6, 3))

for j, y in enumerate(["edge P/A", "edge P/A (reduced)"]):
    ax = axs[j]
    sns.histplot(
        data=df,
        x="block P/A",
        y=y,
        # hue="NZ_JAOSEJ010000001",
        ax=ax,
        legend=False,
    )

sns.despine()
plt.tight_layout()
svfig("7_edges_vs_blocks")
plt.show()

# %%


def df_to_mat(df, val, idx_order):
    sdf = df.pivot(index="si", columns="sj", values=val)
    return sdf.loc[idx_order, idx_order].to_numpy()


M_core = df_to_mat(adf, val="core_div_filtered", idx_order=str_ord)
M_pe = df_to_mat(adf, val="edge_PA", idx_order=str_ord)
M_per = df_to_mat(adf, val="edge_PA_reduced", idx_order=str_ord)
# %%


def matrix_plot(matrix_triplet, svname):
    fig, axs = plt.subplots(
        1,
        4,
        sharey=True,
        figsize=(12, 3.8),
        gridspec_kw={"width_ratios": [0.3, 1, 1, 1]},
    )

    ax = axs[0]
    Phylo.draw(tree, do_show=False, label_func=lambda x: None, axes=ax)

    cmap = plt.get_cmap("viridis_r")

    for i, v in enumerate(matrix_triplet):
        ax = axs[i + 1]
        M, t = v
        g = ax.matshow(M, cmap=cmap)
        ax.set_title(t)
        plt.colorbar(g, ax=ax, shrink=0.8)

    for ax in axs:
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(None)
        ax.set_ylabel(None)

    ax.set_ylim(top=-0.8)
    plt.tight_layout()
    svfig(svname)
    plt.show()


triplet = [
    [M_core, "core-alignment divergence"],
    [M_pe, "edge P/A"],
    [M_per, "edge P/A (reduced)"],
]

matrix_plot(triplet, "tree_vs_edges")

# %%
