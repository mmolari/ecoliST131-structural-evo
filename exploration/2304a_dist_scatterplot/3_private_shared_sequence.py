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
    for suffix in ["png", "pdf", "svg"]:
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
df = adf[mask]

cd = "core genome divergence"
md = "mash distance"
epa = "edge P/A distance"
epar = "edge P/A reduced distance"
bpa = "block P/A distance"
npb = "n. blocks (pariwise projection)"
ps = "private seq. (bp)"
ss = "shared seq. (bp)"
bs = "block_sharing"
es = "edge_sharing"

df = df.rename(
    columns={
        "core_div_filtered": cd,
        "mash_dist": md,
        "edge_PA": epa,
        "edge_PA_reduced": epar,
        "block_PA": bpa,
        "n. blocks": npb,
        "private seq. (bp)": ps,
        "n. shared blocks": ss,
        "n. shared edges": es,
    }
)
# %%

# seaborn pairgrid

g = sns.PairGrid(
    df,
    vars=[cd, md, epa, epar, bpa, npb, ps, ss, bs, es],
    diag_sharey=False,
)
g.map_lower(sns.histplot)
g.map_diag(sns.histplot)
# g.map_upper(sns.kdeplot, fill=True)

# hide upper-diagonal axes
for i, j in zip(*np.triu_indices_from(g.axes, 1)):
    g.axes[i, j].set_visible(False)

plt.tight_layout()
svfig("0_all_comparisons")
plt.show()


# %%

g = sns.PairGrid(
    df,
    vars=[cd, ps, ss],
    diag_sharey=False,
)
g.map_lower(sns.histplot)
g.map_diag(sns.histplot)
# g.map_upper(sns.kdeplot, fill=True)

# hide upper-diagonal axes
for i, j in zip(*np.triu_indices_from(g.axes, 1)):
    g.axes[i, j].set_visible(False)

plt.tight_layout()
svfig("1_priv_shared_seq")
plt.show()

# %%

fig, ax = plt.subplots(1, 1, figsize=(3, 3))
sns.histplot(data=df, x=ps, y=md)
plt.tight_layout()
svfig("2_mash_vs_private_seq")
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

# %%


def matrix_plot(matrix_triplet, svname, revert=None):
    fig, axs = plt.subplots(
        1,
        4,
        sharey=True,
        figsize=(12, 3.8),
        gridspec_kw={"width_ratios": [0.3, 1, 1, 1]},
    )

    if revert is None:
        revert = [False] * len(matrix_triplet)

    ax = axs[0]
    Phylo.draw(tree, do_show=False, label_func=lambda x: None, axes=ax)

    cmap = plt.get_cmap("viridis_r")

    for i, v in enumerate(matrix_triplet):
        ax = axs[i + 1]
        M, t = v
        cm = cmap if not revert[i] else cmap.reversed()
        g = ax.matshow(M, cmap=cm)
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
    [M_ps, "private seq. (bp)"],
    [M_ss, "shared seq. (bp)"],
]

matrix_plot(triplet, "3_coretree_vs_seq", revert=[False, False, True])

# %%

# NZ_JAOSEJ010000001

weird = "NZ_JAOSEJ010000001"
close = "NZ_CP014497"

tree_c = {weird: "C0", close: "C1"}

Phylo.draw(
    tree.root[2][1][0][2][1],
    do_show=False,
    label_func=lambda x: x.name if x.name in str_ord else None,
    label_colors=lambda x: tree_c[x] if x in tree_c else "black",
)
plt.title(None)
sns.despine()
plt.tight_layout()
plt.savefig(f"figs/{weird}--vs--{close}_tree.png", dpi=300, facecolor="white")
plt.show()

# %%


import pypangraph as pp
from pypangraph.pangraph_projector import PanProjector
from pypangraph.visualization_projection import draw_projection
from collections import defaultdict


prefix = "../../results/ST131/pangraph"
pan_file = f"{prefix}/asm20-100-5-polished.json"

pan = pp.Pangraph.load_json(pan_file)
i1, i2 = weird, close

# create projector
ppj = PanProjector(pan)

# project over the pair
pr = ppj.project(i1, i2, exclude_dupl=False)

fig, ax = plt.subplots(1, 1, figsize=(5, 5))

cdict = defaultdict(lambda: plt.get_cmap("rainbow")(np.random.rand()))
draw_projection(
    pr,
    ax=ax,
    color_dict=cdict,
)
ax.legend()
ax.set_xticks([])
ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)

plt.tight_layout()
plt.savefig(f"figs/{weird}--vs--{close}.png", dpi=300, facecolor="white")
plt.show()

# %%
