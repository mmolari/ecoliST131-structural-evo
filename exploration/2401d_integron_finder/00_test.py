# %%
import pandas as pd
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import pathlib

# %%

fld = pathlib.Path("../../results/ST131_ABC/")

fig_fld = pathlib.Path("figs/f00")
fig_fld.mkdir(parents=True, exist_ok=True)

tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]

integron_summary_file = fld / "annotations/integron_finder/integron_summary.tsv"
df = pd.read_csv(integron_summary_file, sep="\t", index_col=0)
df

# %%

cal = df["CALIN"].to_dict()
compl = df["complete"].to_dict()

fig, axs = plt.subplots(
    1, 3, figsize=(8, 10), sharey=True, gridspec_kw={"width_ratios": [1, 0.5, 0.5]}
)

ax = axs[0]


def lab_f(x):
    if x.name not in strains:
        return ""
    if (cal[x.name] > 0) or (compl[x.name] > 0):
        return "~"
    return ""


Phylo.draw(
    tree,
    axes=ax,
    show_confidence=False,
    label_func=lab_f,
    do_show=False,
    label_colors={"~": "red"},
)
ax.set_ylabel("")

ax = axs[1]
ax.barh(np.arange(len(strains)) + 1, [cal[x] for x in strains], height=0.8)
ax.set_xlabel("# CALIN")
ax.set_xticks(np.arange(df["CALIN"].max() + 1))

ax = axs[2]
ax.barh(np.arange(len(strains)) + 1, [compl[x] for x in strains], height=0.8)
ax.set_xlabel("# complete integrons")
ax.set_xticks(np.arange(df["complete"].max() + 1))

for ax in axs:
    ax.grid(axis="y", alpha=0.5)

plt.tight_layout()
plt.savefig(fig_fld / "integron_finder_summary.png")
plt.show()

# %%
