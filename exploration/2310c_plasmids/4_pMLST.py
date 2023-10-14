# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from Bio import Phylo
import json
import argparse
from collections import defaultdict

mlst_df_file = "../../results/ST131_full/plasmids/alleles_summary.csv"
tree_file = "../../results/ST131_full/pangraph/asm20-100-5-filtered-coretree.nwk"
pls_json = "../../config/datasets/ST131_full/plasmids.json"

tree = Phylo.read(tree_file, "newick")
tree.root_at_midpoint()
tree.ladderize()

with open(pls_json, "r") as f:
    pls_dict = json.load(f)

df = pd.read_csv(mlst_df_file, sep=",", index_col=0, dtype=str)
df.fillna("-", inplace=True)
df["typ"] = "F" + df["FII"] + ":A" + df["FIA"] + ":B" + df["FIB"]

# remove untyped
mask = df["typ"] == "F-:A-:B-"
df = df[~mask]

df["typ"].value_counts()
# %%

# P/A matrix for first 8:
I = [l.name for l in tree.get_terminals()]
J = df["typ"].value_counts().index[:8]

pa_df = []
for iso in I:
    res = {"iso": iso}
    if not iso in pls_dict:
        pa_df.append(res)
        continue
    P = pls_dict[iso]
    for plsm in P:
        if not plsm in df.index:
            continue
        t = df.loc[plsm, "typ"]
        if t in J:
            res[t] = True
    pa_df.append(res)
pa_df = pd.DataFrame(pa_df)
pa_df.set_index("iso", inplace=True)
pa_df.fillna(False, inplace=True)
pa_df = pa_df[J]
pa_df = pa_df.loc[I]
pa_df
# %%
# plot tree and matrix
fig, axs = plt.subplots(1, 2, figsize=(10, 10))
ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
ax.set_yticks([])
# ax.set_xticks([])
ax.set_title("ST131 core tree")

ax = axs[1]
sns.heatmap(pa_df, ax=ax, cbar=False, cmap="binary")
ax.set_yticks([])
ax.set_ylabel("Isolate")
ax.set_xlabel("Plasmid type")
ax.set_title("Plasmid types in ST131 isolates")
plt.tight_layout()
sns.despine()
plt.savefig("figs/plasmid_MLST.pdf")
plt.savefig("figs/plasmid_MLST.png", dpi=300, facecolor="white")
plt.show()

# %%
