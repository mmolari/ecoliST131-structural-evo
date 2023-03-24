# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import numpy as np
import json
import pypangraph as pp

from Bio import Phylo
from collections import defaultdict
from pypangraph.pangraph_projector import PanProjector
from pypangraph.visualization_projection import draw_projection

fig_fld = pathlib.Path("figs/pairs")
fig_fld.mkdir(exist_ok=True)


tree = Phylo.read(
    "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk", format="newick"
)
tree.ladderize()

df = pd.read_csv(f"../../results/ST131/distances/summary-asm20-100-5.csv")
df = df.set_index(["si", "sj"])
mask = df.index.get_level_values(0) > df.index.get_level_values(1)
df = df[mask]

info_new = "../../results/ST131/pangraph/asm20-100-5-alignment/filtered_corealignment_info.json"
with open(info_new, "r") as f:
    info = json.load(f)
factor = info["polished aln consensus"] + info["polished aln snps"]

df["corealn snps"] = df["core_div_filtered"] * factor
df["private seq. (kbp)"] = df["private seq. (bp)"] / 1000
df["shared seq. (kbp)"] = df["shared seq. (bp)"] / 2000

# %%

# far
i1 = "NZ_JAOSCJ010000001"
i2 = "NZ_CP076689"
# close
# i1 = "NZ_CP104846"
# i2 = "NZ_CP104848"
# # closish
# i1 = "NZ_SEVU01000007"
# i2 = "NZ_CP035476"
# # intermediate
# i1 = "NZ_CP014497"
# i2 = "NZ_JAOSES010000001"
# # cross-clade
# i1 = "NZ_JAOSCG010000001"
# i2 = "NZ_CP021454"


exc_dupl = False
if i1 < i2:
    i1, i2 = i2, i1

pair = df.loc[(i1, i2)]
pair_idx = df.index.get_loc((i1, i2))

# %%

# load pangenome graph
pangraph_file = "../../results/ST131/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)

# set of strains
strains = pan.strains()

# create projector
ppj = PanProjector(pan)

# project over the pair
pr = ppj.project(i1, i2, exclude_dupl=exc_dupl)

# %%
fig, ax = plt.subplots(1, 1)
draw_projection(
    pr, ax=ax, color_dict=defaultdict(lambda: plt.get_cmap("rainbow")(np.random.rand()))
)
plt.legend()
plt.show()

# %%
# fig, axs = plt.subplots(3, 2, figsize=(8, 8))

fig, axs = plt.subplot_mosaic(
    """
    AADDD
    AADDD
    BBDDD
    BBEEE
    CCEEE
    CCEEE
    """,
    figsize=(8, 8),
)
sns.despine(fig)

ax = axs["A"]
xlab, ylab = "corealn snps", "private seq. (kbp)"
sns.histplot(data=df, x=xlab, y=ylab, ax=ax)
ax.plot([pair[xlab]], [pair[ylab]], "rx")

ax = axs["B"]
xlab, ylab = "corealn snps", "n. blocks"
sns.histplot(data=df, x=xlab, y=ylab, ax=ax)
ax.plot([pair[xlab]], [pair[ylab]], "rx")

ax = axs["C"]
xlab, ylab = "corealn snps", "shared seq. (kbp)"
sns.histplot(data=df, x=xlab, y=ylab, ax=ax)
ax.plot([pair[xlab]], [pair[ylab]], "rx")

ax = axs["D"]
Phylo.draw(
    tree,
    do_show=False,
    axes=ax,
    label_func=lambda x: "" if x.name not in [i1, i2] else x.name,
    label_colors=lambda x: "red",
)
ax.set_ylabel("")
ax.set_yticks([])
ax.spines["left"].set_visible(False)

ax = axs["E"]
draw_projection(
    pr, ax=ax, color_dict=defaultdict(lambda: plt.get_cmap("rainbow")(np.random.rand()))
)
ax.legend()
ax.set_xticks([])
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)


plt.tight_layout()
plt.savefig(fig_fld / f"{i1}-{i2}.pdf")
plt.show()
# %%
