# %%

import os
import pathlib
import json

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import pypangraph as pp
import utils as ut

from Bio import Phylo
from collections import defaultdict, Counter

# %%
"""
plan:
1. create node-named tree
2a clean pahts:
    - remove duplicates
    - remove short blocks (> 500 bp)
    - remove blocks with ambiguous flanking positions
2b plot selected block stats
2c merge transitive edges
3. extract blocks PA
4. extract block context
5. infer normal mugration
6. infer context-aware mugration
"""


data_fld = pathlib.Path("data")
data_fld.mkdir(exist_ok=True)

fig_fld = pathlib.Path("figs/1")
fig_fld.mkdir(exist_ok=True, parents=True)

# %% 1: create node-named tree

tree = ut.load_tree()
for iso, node in enumerate(tree.get_nonterminals()):
    node.name = f"node_{iso:02d}"
Phylo.write(tree, ut.named_nodes_tree_file, "newick", format_branch_length="%.5e")

# %% 2a: clean paths

THR_LEN = 500

pan = ut.load_pangraph()
bdf = pan.to_blockstats_df()
block_Ls = bdf["len"].to_dict()
is_dupl = bdf["duplicated"].to_dict()
is_core = bdf["core"].to_dict()

forb_blocks = pd.read_csv(ut.forb_blocks_file, header=None)[0].to_list()


def keep_block(bid):
    if block_Ls[bid] < THR_LEN:
        return False
    if is_dupl[bid]:
        return False
    if bid in forb_blocks:
        return False
    return True


paths = ut.pangraph_to_nodepath(pan)
paths = ut.clean_up_paths(paths, keep_block)


# %% 2b: plot selected blocks stats

selected_blocks = set()
for p in paths.values():
    selected_blocks.update([n.id for n in p])
selected_blocks = list(selected_blocks)
bdf["category"] = bdf.index.isin(selected_blocks)
bdf["category"] = bdf["category"].map({True: "selected", False: "discarded"})
mask = bdf["core"] & (bdf["category"] == "selected")
bdf.loc[mask, "category"] = "retained core"
# %%

fig, axs = plt.subplots(2, 2, figsize=(10, 8), squeeze=True)


for i, w in enumerate([None, "len"]):
    ax = axs[i, 0]
    bins = np.logspace(1, np.log10(bdf["len"].max()) + 0.1, 30)
    sns.histplot(
        data=bdf,
        x="len",
        hue="category",
        bins=bins,
        weights=w,
        cumulative=False,
        element="step",
        ax=ax,
    )
    ax.set_xscale("log")

    ax = axs[i, 1]
    bins = np.arange(bdf["n. strains"].max() + 1) + 0.5
    sns.histplot(
        data=bdf,
        x="n. strains",
        hue="category",
        bins=bins,
        weights=w,
        cumulative=False,
        element="step",
        ax=ax,
    )

for i in [0, 1]:
    axs[0, i].set_ylabel("n. blocks")
    axs[1, i].set_ylabel("length (bp)")

    axs[i, 0].set_xlabel("block length (bp)")


plt.tight_layout()
plt.savefig(fig_fld / "block_category_distr_filtered.png")
plt.show()

print("n. of blocks per category")
print(bdf["category"].value_counts())
print("weighted by sequence length (Mbp)")
print(bdf.groupby("category")["len"].sum() / 1e6)

# %% 2 compactify paths

mask = [keep_block(bid) for bid in bdf.index]
bdf_sub = bdf[mask].copy().drop(columns=["duplicated", "category", "n. strains"])
pa_sub = pan.to_blockcount_df().T.loc[bdf_sub.index] > 0


def paths_to_edges(paths):
    edges = {}
    for iso, path in paths.items():
        E = []
        L = len(path)
        for l in range(L):
            prev = (l - 1) % L
            edge = ut.Edge(path[prev], path[l])
            E.append(edge)
        edges[iso] = E
    return edges


edges = paths_to_edges(paths)

# all possible edges
E_set = set(sum([E for E in edges.values()], start=[]))

merge_list = {}
for e in E_set:
    pa = {iso: (e in E) for iso, E in edges.items()}
    l, r = e.left.id, e.right.id
    same_l = np.all(pa_sub.loc[l] == pd.Series(pa))
    same_r = np.all(pa_sub.loc[r] == pd.Series(pa))
    if not (same_l and same_r):
        continue
    # print("same", e)
    if (l in merge_list) and (r in merge_list):
        ref = merge_list[l]
        old_ref = merge_list[r]
        for k in merge_list:
            if merge_list[k] == old_ref:
                merge_list[k] = ref
    elif l in merge_list:
        ref = merge_list[l]
    elif r in merge_list:
        ref = merge_list[r]
    else:
        ref = l
    merge_list[l] = ref
    merge_list[r] = ref

# %%

# clean up paths
kept = set(merge_list.values())
deleted = set(merge_list.keys()) - kept
for iso, path in paths.items():
    for l in range(len(path))[::-1]:
        if path[l].id in deleted:
            path.pop(l)

# adapt stats
bdf_sub["parts"] = 1
for k in deleted:
    v = merge_list[k]
    bdf_sub.loc[v, "len"] += bdf_sub.loc[k, "len"]
    bdf_sub.loc[v, "parts"] += 1
    assert bdf_sub.loc[v, "count"] == bdf_sub.loc[k, "count"]
    bdf_sub.drop(k, inplace=True)


# %%
pa_sub.drop(deleted, inplace=True)
# %%

# cleaned up paths
# reduced P/A
# reduced/updated block info
# now redo inference
