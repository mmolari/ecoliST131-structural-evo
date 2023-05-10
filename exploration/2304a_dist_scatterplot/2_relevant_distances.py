# %%
import pathlib
import pandas as pd
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


# tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
# tree = Phylo.read(tree_file, format="newick")
# tree.ladderize()
# str_ord = [x.name for x in tree.get_terminals()]

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

df = df.rename(
    columns={
        "core_div_filtered": cd,
        "mash_dist": md,
        "edge_PA": epa,
        "edge_PA_reduced": epar,
        "block_PA": bpa,
        "n. blocks": npb,
        "private seq. (bp)": ps,
    }
)
# %%


def ax_label(axs):
    for n, ax in enumerate(axs.flatten()):
        ax.text(
            1,
            1,
            chr(65 + n),
            transform=ax.transAxes,
            fontsize="x-large",
            weight="bold",
        )


fig, axs = plt.subplots(2, 2, figsize=(7, 6), sharex=True)
sns.histplot(data=df, x=cd, y=ps, ax=axs[0, 0])
sns.histplot(data=df, x=cd, y=bpa, ax=axs[0, 1])
sns.histplot(data=df, x=cd, y=npb, ax=axs[1, 0])
sns.histplot(data=df, x=cd, y=epa, ax=axs[1, 1])
plt.tight_layout()
sns.despine()
ax_label(axs)
svfig("scatter_1")
plt.show()

# %%
fig, axs = plt.subplots(2, 2, figsize=(6, 6), sharex="col")
sns.histplot(data=df, y=bpa, x=ps, ax=axs[0, 0])
sns.histplot(data=df, y=md, x=ps, ax=axs[1, 0])
sns.histplot(data=df, x=npb, y=epa, ax=axs[0, 1])
sns.histplot(data=df, x=npb, y=epar, ax=axs[1, 1])
plt.tight_layout()
sns.despine()
ax_label(axs)
svfig("scatter_2")
plt.show()

# %%
