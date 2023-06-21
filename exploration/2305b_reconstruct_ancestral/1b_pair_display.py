# %%

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
from Bio import Phylo
import pypangraph as pp

import utils as ut

from collections import defaultdict
from itertools import combinations


# %%

svfld = pathlib.Path("figs/pairs")
svfld.mkdir(exist_ok=True)


# %%

# pair_1 = ["NZ_CP104846", "NZ_CP104848"]
# pair_1 = ["NZ_CP035476", "NZ_CP035720"]
pair_1 = ["NZ_CP076687", "NZ_CP014495"]
# pair_1 = ["NZ_SEVU01000007", "NZ_SEVT01000005"]
pair_2 = ["NZ_CP035516", "NZ_CP035477"]
pairs = pair_1 + pair_2

subfld = svfld / "__".join(pairs)
subfld.mkdir(exist_ok=True)


def svfig(name):
    plt.savefig(subfld / name, dpi=300, facecolor="white", bbox_inches="tight")


# %%

# load tree
tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")

# main_colors = ["magenta", "cyan", "orange", "green"]
main_colors = ["C0", "C1", "C2", "C3"]
iso_colors = {k: main_colors[i] for i, k in enumerate(pairs)}

fig = ut.plot_tree_pairs(tree, pairs, iso_colors)
svfig("pairs_on_tree.png")
plt.show()

# %%
pangraph_file = f"../../results/ST131/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)
# %%

df = pan.to_blockcount_df()

df = df.loc[pairs]

# remove columns which are all zeros
df = df.loc[:, (df != 0).any(axis=0)]

# select columns that are all ones
anchor_genes = df.loc[:, (df == 1).all(axis=0)].columns.to_numpy()

# select duplicated genes: columns that have values >1 at least once
dupl_genes = df.loc[:, (df > 1).any(axis=0)].columns.to_numpy()

# select accessory genes: columns that are neither anchor nor duplicated
acc_genes = df.loc[
    :, ~df.columns.isin(np.concatenate([anchor_genes, dupl_genes]))
].columns.to_numpy()
# %%
ctr = defaultdict(int)
ctr_labels = defaultdict(list)
for idx, row in df[acc_genes].astype(bool).iteritems():
    label = "|".join(sorted(row.index[row]))
    ctr[label] += 1
    ctr_labels[label].append(idx)


# %%

fig = ut.display_overlap(ctr, pairs, iso_colors)
svfig("accessory_overlap.png")
plt.show()

# %%

ctr = defaultdict(int)
ctr_labels = defaultdict(list)
for idx, row in df[dupl_genes].iteritems():
    mask = row > 0
    label = "|".join(sorted(row.index[mask]))
    ctr[label] += 1
    ctr_labels[label].append(idx)

fig = ut.display_overlap(ctr, pairs, iso_colors)
svfig("dupl_overlap.png")
plt.show()

# weird pattern
# %%

# mask = np.all((df > 0).T == np.array([1, 0, 0, 1]), axis=1)
# print(df.T[mask].T)
# q = "SYTUZUAKRF"
# b = pan.blocks[q]
# adf = pan.to_blockcount_df()
# pattern = (adf[q] > 0).to_dict()

# # plot tree, coloring leaves according to pattern
# fig, ax = plt.subplots(figsize=(6, 12))
# Phylo.draw(
#     tree,
#     axes=ax,
#     do_show=False,
#     label_func=lambda x: x.name if x in tree.get_terminals() else "",
#     label_colors=lambda x: "red" if x in pattern and pattern[x] else "black",
#     # branch_labels=lambda x: x.branch_length,
#     # branch_labels_color="red",
#     # branch_labels_size=10,
#     show_confidence=False,
# )
# svfig(f"pattern_{q}.png")
# plt.show()

# %%


block_colors = {k: iso_colors[k] for k in pairs}
for k1, k2 in combinations(pairs, 2):
    lab = "|".join(sorted([k1, k2]))
    block_colors[lab] = np.mean(
        [mpl.colors.to_rgb(iso_colors[k]) for k in [k1, k2]], axis=0
    )
for k1, k2, k3 in combinations(pairs, 3):
    lab = "|".join(sorted([k1, k2, k3]))
    block_colors[lab] = np.mean(
        [mpl.colors.to_rgb(iso_colors[k]) for k in [k1, k2, k3]],
        axis=0,
    )
block_colors["|".join(sorted(pairs))] = "lightgrey"


fig, ax = plt.subplots(1, 1, figsize=(20, 20), subplot_kw={"projection": "polar"})

bl_lens = pan.to_blockstats_df()["len"].to_dict()
core_anchor = anchor_genes[np.argmax([bl_lens[x] for x in anchor_genes])]

Lmax = max([np.sum([bl_lens[b] for b in pan.paths[k].block_ids]) for k in pairs])
arc = 2 * np.pi / Lmax

for i, k in enumerate(pairs):
    y = 5 + (i * 0.25)
    p = pan.paths[k]
    bls = list(p.block_ids)
    anch = bls.index(core_anchor)
    bls = bls[anch:] + bls[:anch]
    if not p.block_strands[anch]:
        bls = bls[::-1]
    l = 0
    for b in bls:
        x = np.array([l, l + bl_lens[b]]) * arc
        strains = df.index[df[b] > 0]
        dy = len(strains) * 0.05
        l += bl_lens[b]
        lab = "|".join(sorted(strains))
        c = block_colors[lab]
        ax.plot(x, [y + dy, y + dy], color=c)

ax.grid(False)
ax.set_ylim(bottom=0)
ax.set_xticks([])
ax.set_yticks([])
svfig("polar_full.png")
plt.show()
