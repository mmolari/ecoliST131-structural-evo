# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import Phylo
import pathlib

# %%
res_fld = pathlib.Path("../../results/ST131_ABC/annotations/defense_finder")
df_file = res_fld / "systems_summary.tsv"
Sdf = pd.read_csv(df_file, sep="\t")
Sdf["iso"] = Sdf["sys_beg"].str.split("_").apply(lambda x: "_".join(x[:-1]))

df_file = res_fld / "genes_summary.tsv"
Gdf = pd.read_csv(df_file, sep="\t")

tree_file = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [l.name for l in tree.get_terminals()]

fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(exist_ok=True, parents=True)

# %%
fig, axs = plt.subplots(1, 2, figsize=(9, 4))
ax = axs[0]
sns.histplot(Sdf["iso"].value_counts(), ax=ax, discrete=True)
ax.set_xlabel("n. systems")
ax.set_ylabel("n. isolates")

ax = axs[1]
sns.histplot(Sdf["genes_count"], ax=ax, discrete=True)
ax.set_xlabel("n. genes")
ax.set_ylabel("n. systems")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "system_histograms.png")
plt.show()
# %%

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
vc = Sdf["subtype"].value_counts().iloc[:20]
sns.barplot(vc)
ax.set_xticklabels(ax.get_xticklabels(), rotation=40, horizontalalignment="right")
plt.tight_layout()
sns.despine()
plt.savefig(fig_fld / "system_types_hist.png")
plt.show()

# %%

fig, axs = plt.subplots(
    1, 3, figsize=(10, 10), sharey=True, gridspec_kw={"width_ratios": [1, 4, 1]}
)

Ys = np.arange(len(strains)) + 1

ax = axs[0]
Phylo.draw(tree, do_show=False, show_confidence=False, axes=ax, label_func=lambda x: "")

ax = axs[1]
vc = Sdf["subtype"].value_counts().iloc[:18]
i = 0
for k, v in vc.items():
    mask = Sdf["subtype"] == k
    iso_ct = Sdf[mask]["iso"].value_counts().to_dict()
    Ns = np.array([iso_ct[s] if s in iso_ct else 0 for s in strains])
    print(k, Ns.max())
    delta = Ns / 5
    ax.barh(y=Ys, width=delta, left=i, height=1.0, color="gray")
    i += 1

ax.set_xticks(np.arange(len(vc)))
ax.set_xticklabels(vc.index, rotation=90, horizontalalignment="center")
ax.grid(axis="x", alpha=0.4)


ax = axs[2]
ns = Sdf["iso"].value_counts()
ax.barh(y=Ys, width=ns.loc[strains].to_numpy(), height=1.0, color="royalblue")
ax.set_yticks(np.arange(0, len(strains), 25)[1:])
ax.set_xlabel("tot n. systems")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "systems_tree_distr.png")
plt.savefig(fig_fld / "systems_tree_distr.pdf")
plt.show()

# %%
