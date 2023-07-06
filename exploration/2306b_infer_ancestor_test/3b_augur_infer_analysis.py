# %%

import pathlib
import os
import json

import itertools as itt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import pypangraph as pp

from Bio import Phylo
from collections import Counter

fld = pathlib.Path("data/mugration")
fld.mkdir(exist_ok=True)

fig_fld = pathlib.Path("figs/mugration")
fig_fld.mkdir(exist_ok=True)

# %%

prefix = "../../results/ST131/pangraph"
original_tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
tree_file = fld / "tree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

# load tree and pangraph
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()

# %%

# select non-duplicated accessory genes
mask = ~bdf["core"]
mask &= ~bdf["duplicated"]

# list of accessory blocks ids
B_acc = bdf[mask].index.to_list()

pa_inference_file = fld / "infer_pa.json"
with open(pa_inference_file, "r") as f:
    pa_inference = json.load(f)


# %%

summary = []
for b in B_acc:
    info = pa_inference[b]
    pa = info["pa_pattern"]["pa"]
    r = info["rates"]
    ev = info["pa_pattern"]["events"]
    # count event type
    n_gain = len([e for e in ev if e[1] == "gain"])
    n_loss = len([e for e in ev if e[1] == "loss"])
    # count n. of leaves
    pa_count = [v for k, v in pa.items() if not k.startswith("node_")]
    pa_count = Counter(pa_count)
    summary.append(
        {
            "bid": b,
            "n_events": len(ev),
            "n_gains": n_gain,
            "n_losses": n_loss,
            "n_leaves_P": pa_count["P"],
            "n_leaves_A": pa_count["A"],
        }
    )
df = pd.DataFrame(summary).set_index("bid")
df["max_events"] = df[["n_leaves_P", "n_leaves_A"]].min(axis=1)
Ls = bdf["len"].to_dict()
df["len"] = df.apply(lambda x: Ls[x.name], axis=1)
df

# %%
df[["n_events", "n_gains", "n_losses", "max_events"]].value_counts()

# %%

fig, axs = plt.subplots(1, 2, figsize=(10, 4))

ax = axs[0]
M = max(df["n_events"].max(), df["max_events"].max())
bins = np.arange(M + 2) - 0.5
bins = (bins, bins)
sns.histplot(
    df,
    x="n_events",
    y="max_events",
    bins=bins,
    cbar=True,
    cbar_kws={
        # "shrink": 0.5,
        "label": "n. blocks",
    },
    norm=mpl.colors.LogNorm(),
    vmin=None,
    vmax=None,
    # weights="len",
    ax=ax,
)
ax.set_xlabel("n. events")
ax.set_ylabel("n. leaves with minority pattern")
ax.plot([0, M], [0, M], "k--")

ax = axs[1]
M = max(df["n_gains"].max(), df["n_losses"].max())
bins = np.arange(M + 2) - 0.5
bins = (bins, bins)
sns.histplot(
    df,
    x="n_gains",
    y="n_losses",
    bins=bins,
    cbar=True,
    cbar_kws={
        # "shrink": 0.5,
        "label": "n. blocks",
    },
    norm=mpl.colors.LogNorm(),
    vmin=None,
    vmax=None,
    # weights="len",
    ax=ax,
)
ax.set_xlabel("n. gains")
ax.set_ylabel("n. losses")
ax.plot([0, M], [0, M], "k--")
plt.tight_layout()
plt.savefig(fig_fld / "n_events.pdf")
plt.show()
# %%

mask = df["max_events"] < df["n_events"]
df[mask]

# %%
bid = "WNAPVOJCVV"
# bid = "JGMLMOYQKZ"
bid = "MANREFTSIL"

tree = Phylo.read(tree_file, "newick")

info = pa_inference[bid]
pa = info["pa_pattern"]["pa"]
ev = info["pa_pattern"]["events"]


def color_tree(node):
    # if node.is_terminal():
    #     if pa[node.name] == "P":
    #         node.color = "green"
    #     else:
    #         node.color = "red"
    # else:
    for nn, et in ev:
        if node.name == nn:
            node.color = "lime" if et == "gain" else "red"
    if node.color is None:
        node.color = "black"
    for c in node.clades:
        color_tree(c)


def label_tree(node):
    if node.is_terminal():
        return node.name
    else:
        return ""


def lab_colors(nn):
    if len(nn) == 0:
        return None
    if pa[nn] == "P":
        return "green"
    else:
        return "white"


color_tree(tree.root)

fig, ax = plt.subplots(1, 1, figsize=(10, 12))
Phylo.draw(tree, label_func=label_tree, label_colors=lab_colors, axes=ax, do_show=False)
plt.tight_layout()
plt.show()


# %%
# tyrosine-type recombinase/integrase, ~ 1kb, almost exact overlap
pan.blocks["WNAPVOJCVV"].sequence
# %%

# DNA transfer protein
pan.blocks["JGMLMOYQKZ"].sequence
# %%

# TODO:
# - document this
# - look for branches with good support (1 event / non-trivial)
# - what is the rest?
# - hist of accessory blocks PA frequencies

# %%
sdf = bdf.loc[B_acc]["n. strains"]
sns.histplot(sdf, bins=np.arange(sdf.max() + 2) - 0.5)
plt.show()
# %%
df.loc["MANREFTSIL"]
# %%
