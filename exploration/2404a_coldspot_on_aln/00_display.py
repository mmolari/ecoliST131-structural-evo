# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo
import pathlib

fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(fname, "newick")
strains = [x.name for x in tree.get_terminals()]

fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/j_aln_pos/aln_pos.csv"
aln_pos = pd.read_csv(fname, index_col=0)
fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/j_aln_pos/junct_pos.csv"
junct_pos = pd.read_csv(fname, index_col=0)
fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/j_aln_pos/SNPs.csv"
snps_info = pd.read_csv(fname, index_col=0)

fname = "../../results/ST131_ABC/rates/asm20-100-5/terminal_coldspot.csv"
csdf = pd.read_csv(fname, index_col=0)

# %%
fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
js = pd.read_csv(fname, index_col=0)
js = js[js.index.isin(junct_pos.index)]
js.index.name = "junction"
js = js.join(junct_pos, on="junction", validate="one_to_one")

fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_pangenome.csv"
ep = pd.read_csv(fname, index_col=0)
ep = ep[ep.index.isin(js.index)]
ep.index.name = "junction"
js = js.join(ep, on="junction", validate="one_to_one")


# %%
fig, axs = plt.subplots(3, 1, figsize=(20, 10), sharex=True)
ax = axs[0]
sns.histplot(
    snps_info.index,
    binwidth=10000,
    element="step",
    ax=ax,
)

ax = axs[1]
sns.scatterplot(
    js,
    x="position",
    y="n_categories",
    hue=js["n_iso"] < len(strains),
    alpha=0.5,
    ax=ax,
)
# legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=["conserved", "broken"], title="synteny")
ax.set_yscale("log")
ax.set_ylabel("n. path categories")

ax = axs[2]
sns.histplot(
    js,
    x="position",
    weights="pangenome_len",
    binwidth=10000,
    element="step",
    ax=ax,
)
# ax.set_yscale("log")
ax.set_ylabel("pangenome length")

for ax in axs:
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(fig_fld / "snps_and_junctions.png")
plt.show()


# %%

# most frequent entry in SNPs matrix:
consensus = snps_info.apply(lambda x: x.value_counts().idxmax(), axis=1)
is_consensus = snps_info.T == consensus
# %%

fig, axs = plt.subplots(
    1, 2, figsize=(20, 10), sharey=True, gridspec_kw={"width_ratios": [0.1, 1]}
)

ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")

ax = axs[1]
for n, strain in enumerate(strains):
    ic = is_consensus.loc[strain]
    snps = ic.index[np.argwhere(~ic).flatten()]
    ax.scatter(snps, np.full_like(snps, n + 1), marker=".", alpha=0.1, color="k")

for j, row in js.iterrows():
    pos = row["position"]

    if row["n_iso"] < len(strains):
        # ax.axvline(pos, color="r", alpha=0.5)
        continue

    if row["singleton"]:
        iso = csdf.loc[j]["event_iso"]
        iso_y = strains.index(iso) + 1
        tp = csdf.loc[j]["event_type"]

        color = {
            "gain": "g",
            "loss": "r",
            "other": "y",
        }

        c = color.get(tp, "b")

        # def get_color(row):
        #     if row["genomad"] > 0:
        #         return "C0"
        #     if row["defensefinder"] > 0:
        #         return "C3"
        #     if row["isescan"] > 0:
        #         return "C1"

        # c = get_color(csdf.loc[j])

        ax.plot([pos], [iso_y], color=c, marker="x")


plt.tight_layout()
plt.savefig(fig_fld / "aln_matrix_and_coldspots.png")
plt.show()


# %%
