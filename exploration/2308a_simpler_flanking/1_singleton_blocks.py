# %%

# looking for single accessory blocks always flanked by core blocks

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import utils as ut

from collections import defaultdict
from Bio import Phylo

svfld = ut.fig_fld / "simple_flanking"
svfld.mkdir(exist_ok=True, parents=True)

# %%

pan = ut.load_pangraph()
paths = ut.pangraph_to_path_dict(pan)
bdf = pan.to_blockstats_df()
is_core = bdf["core"]

# %%
Js = {iso: ut.to_junctions(path.nodes) for iso, path in paths.items()}

# %%

# look for non-duplicated accessory blocks
mask = (~bdf["core"]) & (~bdf["duplicated"])
candidates = bdf[mask].index.to_numpy()
print(f"n. accessory candidates = {len(candidates)}")

# %%

# check if they are flanked by core blocks
accepted = list(candidates)

for iso, J in Js.items():
    for j in J:
        bid = j.center.id
        if not bid in accepted:
            continue
        l, r = j.left.id, j.right.id
        if not is_core[l] or not is_core[r]:
            accepted.remove(bid)

# %%

sns.scatterplot(
    data=bdf.loc[accepted],
    x="len",
    y="count",
)
plt.xscale("log")
plt.xlabel("block length (bp)")
plt.ylabel("block count")
plt.tight_layout()
plt.savefig(svfld / "accepted_block_len_vs_count.png")
plt.show()

# %%

# check how many contexts
contexts = defaultdict(dict)
for iso, J in Js.items():
    for j in J:
        bid = j.center.id
        if not bid in accepted:
            continue
        contexts[bid][iso] = j

for bid, C in contexts.items():
    print(f"bid = {bid}, n. contexts = {len(set(C.values()))}")

# %%

# presence/absence pattern
PA = pan.to_blockcount_df()
PA = PA[accepted]

mask = (PA.sum(axis=0) > 1) & (PA.sum(axis=0) < (len(pan.strains()) - 1))
non_singleton = PA.columns[mask]
# %%

for bid in non_singleton:
    pa = PA[bid]
    tree = ut.load_tree()
    fig, ax = plt.subplots(figsize=(10, 10))
    Phylo.draw(
        tree,
        label_func=lambda x: x.name if x in tree.get_terminals() else None,
        label_colors=lambda x: "C2" if pa[x] else "C3",
        do_show=False,
        show_confidence=False,
        axes=ax,
    )
    plt.title(f"bid = {bid}, L = {bdf.loc[bid, 'len']} bp")
    plt.tight_layout()
    plt.savefig(svfld / f"tree_{bid}.png", facecolor="white")
    plt.show()

# %%
