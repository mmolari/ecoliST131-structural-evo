# %%
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from Bio import Phylo
import numpy as np
from collections import defaultdict

# %%
tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
allele_files = {
    "fimH": "../../results/ST131/assembly_qc/alleles/summary/fimH.tsv",
    "fimH_eb": "../../results/ST131/assembly_qc/alleles/summary/fimH_eb.tsv",
    "gyrA_eb": "../../results/ST131/assembly_qc/alleles/summary/gyrA_eb.tsv",
    "parC_eb": "../../results/ST131/assembly_qc/alleles/summary/parC_eb.tsv",
}

# %%
tree = Phylo.read(tree_file, "newick")
leaves = [l.name for l in tree.get_terminals()]
y_loc = {l: i + 1 for i, l in enumerate(leaves)}

dfs = {locus: pd.read_csv(f, sep="\t") for locus, f in allele_files.items()}

# %%


# color generator
def cgen():
    cm = plt.cm.get_cmap("tab20")
    for i in range(20):
        yield cm(i)


def colordict():
    cg = cgen()
    cd = defaultdict(lambda: cg.__next__())
    return cd


fig, axs = plt.subplots(
    1, 3, figsize=(5, 10), sharey=True, gridspec_kw={"width_ratios": [3, 1, 1]}
)

ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")

ax = axs[1]
x = 0
xticks = []
xticklabels = []
cdict = colordict()
for locus, df in dfs.items():
    df = df.set_index("iso")
    for iso in leaves:
        match = df.loc[iso, "match"]
        if match != "+":
            continue
        y = y_loc[iso]
        allele = int(df.loc[iso, "allele"])
        allele = f"{locus}_{allele}"
        c = cdict[allele]
        ax.scatter(x, y, color=c, marker="s", edgecolor="black")
    xticks.append(x)
    xticklabels.append(locus)
    x += 1
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, rotation=90)
ax.spines["left"].set_visible(False)
ax.set_xlim(-0.5, x - 0.5)

for ax in axs:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.5)

ax = axs[2]
y = 2
for key, val in cdict.items():
    ax.text(0, y, key, color=val)
    y += 2
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.grid(False)
ax.set_xticks([])

plt.tight_layout()
plt.savefig("figs/alleles.png", dpi=300, facecolor="white")
plt.show()

# %%
