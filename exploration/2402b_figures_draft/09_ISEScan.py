# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
import pathlib

fig_fld = pathlib.Path("figs/f09")
fig_fld.mkdir(parents=True, exist_ok=True)

tree_file = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]

# %%
df_file = "../../results/ST131_ABC/annotations/isescan/is_summary.tsv"
df = pd.read_csv(df_file, sep="\t")
df["type"] = df["type"].map({"c": "complete", "p": "partial"})
df.head()

fam_order = df["family"].value_counts().index.to_list()
df["family"] = pd.Categorical(df["family"], categories=fam_order)
# %%

fig, axs = plt.subplots(1, 2, figsize=(8, 4), gridspec_kw={"width_ratios": [3, 2]})
ax = axs[0]
sns.histplot(df, x="family", ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_xlabel("IS Family")
ax.set_ylabel("count")

ax = axs[1]
sns.histplot(df["seqID"].value_counts())
ax.set_xlabel("n. of ISs")
ax.set_ylabel("n. of isolates")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "is_hist.pdf")
plt.show()

# %%

sdf = df.groupby(["family", "seqID"], observed=False).size().reset_index(name="counts")
sdf = sdf.pivot(index="seqID", columns="family", values="counts")


fig, axs = plt.subplots(
    1, 3, figsize=(9, 10), sharey=True, gridspec_kw={"width_ratios": [1, 4, 1]}
)

ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")

ax = axs[1]
Y = sdf.max().max()
for nf, fam in enumerate(sdf.columns):
    y = sdf[fam][strains].to_numpy()
    ax.barh(
        np.arange(len(y)) + 1,
        # np.log(y + 1, where=y > 0, out=0.0 * y) / 4,
        # np.sqrt(y) / np.sqrt(Y),
        y / Y,
        color="gray",
        height=1,
        left=nf,
    )


ax.set_yticks([])
ax.set_xlim(0, len(sdf.columns))
ax.set_xticks(np.arange(len(sdf.columns)))
ax.set_xticklabels(sdf.columns, rotation=90)
ax.grid(True, axis="x", alpha=0.2)
ax.set_xlabel("IS Family")

ax = axs[2]
plt.barh(
    np.arange(len(strains)) + 1,
    df["seqID"].value_counts()[strains],
    color="gray",
    height=1,
)

plt.tight_layout()
plt.savefig(fig_fld / "is_tree.pdf")
plt.show()


# %%

# %%
