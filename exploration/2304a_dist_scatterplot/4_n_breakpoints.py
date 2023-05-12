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
    "n. pairwise blocks",
    "block P/A",
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
svfig("4_n_blocks")
plt.show()

# %%

# three seaborn histograms

fig, axs = plt.subplots(1, 3, figsize=(9, 3))


sns.histplot(
    data=df,
    x="n. pairwise blocks",
    y="n. breakpoints",
    hue="NZ_JAOSEJ010000001",
    ax=axs[0],
    # legend=False,
)
sns.move_legend(axs[0], "lower right", bbox_to_anchor=(1.1, 0.0), ncol=2)


sns.histplot(
    data=df,
    x="private seq. (bp)",
    y="block P/A",
    hue="NZ_JAOSEJ010000001",
    ax=axs[1],
    legend=False,
)
# sns.move_legend(axs[1], "lower right", bbox_to_anchor=(1.1, 0.0), ncol=2)


sns.histplot(
    data=df,
    x="private seq. (bp)",
    y="n. pairwise blocks",
    hue="NZ_JAOSEJ010000001",
    ax=axs[2],
    legend=False,
)
# sns.move_legend(axs[2], "lower right", bbox_to_anchor=(1.1, 0.0), ncol=2)

sns.despine()
plt.tight_layout()
svfig("5_blocks_vs_priv_seq")
plt.show()

# %%
