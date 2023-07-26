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
from collections import defaultdict

# %%
"""
plan:
1. create node-named tree
2a clean pahts:
    - remove duplicates
    - remove short blocks (> 500 bp)
2b plot selected block stats
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

# clean up paths
paths = {}
for p in pan.paths:
    B = p.block_ids
    S = p.block_strands
    path = []
    for b, s in zip(B, S):
        if block_Ls[b] < THR_LEN:
            continue
        if is_dupl[b]:
            continue
        node = ut.Node(b, s)
        path.append(node)
    paths[p.name] = path

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
plt.savefig(fig_fld / "block_category_distr.png")
plt.show()

print("n. of blocks per category")
print(bdf["category"].value_counts())
print("weighted by sequence length (Mbp)")
print(bdf.groupby("category")["len"].sum() / 1e6)

# %% 3. extract block PA
PA = pan.to_blockcount_df()
PA = PA.loc[:, selected_blocks] > 0
for bid, col in PA.iteritems():
    PA.loc[:, bid] = PA.loc[:, bid].map({True: "P", False: "-"})
PA.index.name = "#name"
PA.to_csv(data_fld / "PA_simple.csv")

# %% 4. extract block context

path_adj = {}
for iso, p in paths.items():
    path_adj[iso] = ut.to_core_adjacencies(p, is_core)
path_adj


PA_ctx = PA.copy()
PA_ctx.loc[:, :] = "-"

adj_letter = defaultdict(lambda: {"next": "A"})
for iso, p in path_adj.items():
    for b, j in p.items():
        if not (j in adj_letter[b]):
            adj_letter[b][j] = adj_letter[b]["next"]
            next_letter = chr(ord(adj_letter[b]["next"]) + 1)
            adj_letter[b]["next"] = next_letter
        PA_ctx.loc[iso, b] = adj_letter[b][j]

PA_ctx.index.name = "#name"
PA_ctx.to_csv(data_fld / "PA_context.csv")

# %% 5: infer normal mugration

mugr_fld = data_fld / "mugration"
mugr_fld.mkdir(exist_ok=True)


pa_inference = {}

mask = bdf["category"] == "selected"
Bs = bdf[mask].index.to_list()
for i, bid in enumerate(Bs):
    print(f"processing {bid} \t - \t {i+1}/{len(Bs)}")
    out_dir = mugr_fld / bid

    # run treetime mugration
    states = data_fld / "PA_simple.csv"

    cmd = f"""
    treetime mugration \
        --tree {ut.named_nodes_tree_file} \
        --states {states} \
        --attribute {bid} \
        --outdir {out_dir} \
    """

    os.system(cmd)

    # output file names
    nexus_file = out_dir / "annotated_tree.nexus"
    rates_file = out_dir / "GTR.txt"

    # parse rates
    rates = ut.parse_gtr(out_dir / "GTR.txt")

    # parse inferred presence/absence pattern
    pa_pattern = ut.parse_nexus(nexus_file)

    os.system(f"rm -r {out_dir}")

    pa_inference[bid] = {
        "pa_pattern": pa_pattern,
        "rates": rates,
    }

# write to json
pa_inference_file = data_fld / "infer_pa_simple.json"
with open(pa_inference_file, "w") as f:
    json.dump(pa_inference, f)

# %%


pa_inference = {}

mask = bdf["category"] == "selected"
Bs = bdf[mask].index.to_list()
for i, bid in enumerate(Bs):
    print(f"processing {bid} \t - \t {i+1}/{len(Bs)}")
    out_dir = mugr_fld / bid

    # run treetime mugration
    states = data_fld / "PA_context.csv"

    cmd = f"""
    treetime mugration \
        --tree {ut.named_nodes_tree_file} \
        --states {states} \
        --attribute {bid} \
        --outdir {out_dir} \
    """

    os.system(cmd)

    # output file names
    nexus_file = out_dir / "annotated_tree.nexus"
    rates_file = out_dir / "GTR.txt"

    # parse rates
    rates = ut.parse_gtr(out_dir / "GTR.txt")

    # parse inferred presence/absence pattern
    pa_pattern = ut.parse_nexus(nexus_file)

    os.system(f"rm -r {out_dir}")

    pa_inference[bid] = {
        "pa_pattern": pa_pattern,
        "rates": rates,
    }

# write to json
pa_inference_file = data_fld / "infer_pa_context.json"
with open(pa_inference_file, "w") as f:
    json.dump(pa_inference, f)

# %%
