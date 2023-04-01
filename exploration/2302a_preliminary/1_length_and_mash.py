# %%
import pathlib
import numpy as np
import pandas as pd
from Bio import SeqIO, Phylo
import matplotlib.pyplot as plt
import seaborn as sns

fig_fld = pathlib.Path("figs")
fig_fld.mkdir(exist_ok=True)


def svfig(svname):
    for suffix in ["png", "pdf", "svg"]:
        plt.savefig(fig_fld / f"{svname}.{suffix}", dpi=300)


# %%
df = pd.read_csv("../../results/ST131/distances/mash_dist.csv")
mash_df = df.pivot(index="si", columns="sj", values="mash_dist")

tree = Phylo.read(
    "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk", format="newick"
)
tree.ladderize()
order = [x.name for x in tree.get_terminals()]
N = len(order)

mash_df = mash_df[order].loc[order]

# %%
fig, axs = plt.subplots(
    1, 2, figsize=(10, 6), gridspec_kw={"width_ratios": [1, 5]}, sharey=True
)

ax = axs[0]
Phylo.draw(tree, do_show=False, axes=ax, label_func=lambda x: None)
for k in ["top", "right", "left"]:
    ax.spines[k].set_visible(False)
ax.set_xlabel("")
ax.set_ylabel("")

ax = axs[1]
g = ax.matshow(mash_df)
plt.colorbar(g, ax=ax, shrink=0.8, label="mash distance")

ax.set_ylim(top=-1)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout()
svfig("mash_dist")
plt.show()
# %%
Ls = {}
for o in order:
    s = SeqIO.read(f"../../data/fa/{o}.fa", format="fasta")
    Ls[o] = len(s) / 1e6
Ls = pd.Series(Ls, name="len")


plt.scatter([Ls[o] for o in order], np.arange(N))
plt.show()

print(f"mean : {Ls.mean()}")
print(f"std : {Ls.std()}")

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 1.7))
sns.swarmplot(data=Ls, orient="h", ax=ax)
for k in ["top", "right", "left"]:
    ax.spines[k].set_visible(False)
# ax.axvline(Ls["NZ_JAOSEJ010000001"])
ax.set_yticks([])
ax.set_xlabel("genome length (Mbp)")
plt.tight_layout()
svfig("len_distr")
plt.show()
# %%
