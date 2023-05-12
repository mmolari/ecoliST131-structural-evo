# %%
import pathlib
import pandas as pd
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

# %%
df_file = "../../results/ST131/distances/summary-asm20-100-5.csv"
adf = pd.read_csv(df_file)
mask = adf["si"] > adf["sj"]
df = adf[mask]
# %%
cols = [
    "si",
    "sj",
    "core_div_naive",
    "mash_dist",
    "private seq. (bp)",
    "shared seq. (bp)",
    "n. breakpoints",
    "part. entropy",
    "n. blocks",
    "core_div_filtered",
    "edge_PA",
    "edge_sharing",
    "block_PA",
    "block_sharing",
]
# %%
sns.histplot(data=df, x="core_div_filtered", y="shared seq. (bp)")
plt.show()
sns.histplot(data=df, x="core_div_filtered", y="private seq. (bp)")
plt.show()
sns.histplot(data=df, x="core_div_filtered", y="edge_sharing")
plt.show()
sns.histplot(data=df, x="core_div_filtered", y="edge_PA")
plt.show()
sns.histplot(data=df, x="core_div_filtered", y="edge_PA_reduced")
plt.show()
sns.histplot(data=df, x="n. blocks", y="edge_PA")
plt.show()
sns.histplot(data=df, x="n. blocks", y="edge_PA_reduced")
plt.show()
# %%


def df_to_mat(df, val, idx_order):
    sdf = df.pivot(index="si", columns="sj", values=val)
    return sdf.loc[idx_order, idx_order].to_numpy()


M_core = df_to_mat(adf, val="core_div_filtered", idx_order=str_ord)
M_ss = df_to_mat(adf, val="shared seq. (bp)", idx_order=str_ord)
M_se = df_to_mat(adf, val="edge_sharing", idx_order=str_ord)
M_ps = df_to_mat(adf, val="private seq. (bp)", idx_order=str_ord)
M_pe = df_to_mat(adf, val="edge_PA", idx_order=str_ord)
M_per = df_to_mat(adf, val="edge_PA_reduced", idx_order=str_ord)
M_sb = df_to_mat(adf, val="block_sharing", idx_order=str_ord)
M_pb = df_to_mat(adf, val="block_PA", idx_order=str_ord)
M_pwb = df_to_mat(adf, val="n. blocks", idx_order=str_ord)
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
    [M_ss, "shared seq. (bp)"],
    [M_se, "shared edges (n)"],
]

matrix_plot(triplet, "tree_vs_shared")

# %%
triplet = [
    [M_core, "core-alignment divergence"],
    [M_ps, "private seq. (bp)"],
    [M_pe, "edge P/A dist (n)"],
]
matrix_plot(triplet, "tree_vs_private")

# %%

triplet = [
    [M_core, "core-alignment divergence"],
    [M_pb, "block P/A dist"],
    [M_pe, "edge P/A dist (n)"],
]
matrix_plot(triplet, "block_vs_edge_PA")


triplet = [
    [M_core, "core-alignment divergence"],
    [M_sb, "shared blocks (n)"],
    [M_se, "shared edges (n)"],
]
matrix_plot(triplet, "block_vs_edge_shared")

# %%

triplet = [
    [M_core, "core-alignment divergence"],
    [M_pb, "block P/A"],
    [M_pwb, "n. pairwise blocks"],
]
matrix_plot(triplet, "blocks_vs_core_tree")

# %%

triplet = [
    [M_core, "core-alignment divergence"],
    [M_pe, "edge P/A dist (n)"],
    [M_per, "edge P/A dist reduced (n)"],
]
matrix_plot(triplet, "tree_vs_edges")

# %%
