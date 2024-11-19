# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pathlib
from collections import defaultdict
import pypangraph as pp
import itertools as itt

# look into inversion that happened multiple times
# 115 kbps
# isolates:
# NZ_CP107114.1
# NZ_CP107184.1
# NZ_CP107182.1

fld = pathlib.Path("../../results/ST131_ABC")
fig_fld = pathlib.Path("figs/f00b")
fig_fld.mkdir(exist_ok=True, parents=True)

colors_file = fld / "pangraph/coresynt-asm20-100-5/blocks.csv"
mergers_file = fld / "pangraph/coresynt-asm20-100-5/mergers.csv"

cl_df = pd.read_csv(colors_file, index_col=0)
mg_df = pd.read_csv(mergers_file, index_col=0)

pangraph_file = fld / "pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
# %%

positions = [
    {"iso": "NZ_CP107114.1", "color": "C0", "start": 2.45e6, "end": 3.05e6},
    {"iso": "NZ_CP107184.1", "color": "C1", "start": 4.3e6, "end": 4.9e6},
    {"iso": "NZ_CP107182.1", "color": "C2", "start": 0.2e6, "end": 0.8e6},
    {"iso": "NZ_CP032201.1", "color": "C3", "start": 0e6, "end": 0.6e6},
    {"iso": "NZ_LS992180.1", "color": "C4", "start": 4.2e6, "end": 4.8e6},
]
positions = pd.DataFrame(positions).set_index("iso")


def position_dict(path):
    pos = path.block_positions
    blocks = path.block_ids
    strand = path.block_strands
    pos_dict = defaultdict(list)
    for i, b in enumerate(blocks):
        pos_dict[b].append((pos[i], pos[i + 1], strand[i]))
    return pos_dict


def in_interval(i, lims):
    S, E = lims
    return (i >= S) and (i <= E)


def dotplot(pan, isoA, isoB, bdf, cl_df, mg_df, ax, xlims, ylims):
    pA = pan.paths[isoA]
    pB = pan.paths[isoB]

    posA = position_dict(pA)
    posB = position_dict(pB)

    for b, PA in posA.items():
        if b not in posB:
            continue
        for bA, eA, sA in PA:
            for bB, eB, sB in posB[b]:
                if (bA - eA) > 2e6 or (bB - eB) > 2e6:
                    continue
                if not (in_interval(bA, xlims) or in_interval(eA, xlims)):
                    continue
                if not (in_interval(bB, ylims) or in_interval(eB, ylims)):
                    continue
                x = [bA, eA]
                y = [bB, eB] if sA == sB else [eB, bB]
                col = "gray"
                lw = 0.5
                if bdf.loc[b, "core"] and bdf.loc[b, "len"] > 500:
                    if b in mg_df.index:
                        mg = mg_df.loc[b, "0"]
                    else:
                        mg = b
                    col = cl_df.loc[mg, "color"]
                    lw = 3
                    # lw = 1
                    # if mg == "RYYAQMEJGY":
                    #     lw = 3
                ax.plot(x, y, color=col, lw=lw)
    ax.set_aspect("equal")


N_iso = len(positions)
fig, axs = plt.subplots(
    N_iso - 1,
    N_iso - 1,
    figsize=(13, 13),
    sharex="col",
    sharey="row",
)

for i, j in itt.product(range(N_iso), range(N_iso)):
    ax_x = i
    ax_y = j - 1
    if i >= j:
        if (i < N_iso - 1) and (j > 0):
            ax = axs[ax_y, ax_x]
            ax.set_visible(False)
        continue

    ax = axs[ax_y, ax_x]
    iso_x = positions.index[i]
    iso_y = positions.index[j]

    xlims = positions.loc[iso_x][["start", "end"]].to_list()
    ylims = positions.loc[iso_y][["start", "end"]].to_list()
    dotplot(pan, iso_x, iso_y, bdf, cl_df, mg_df, ax, xlims, ylims)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    if ax_x == 0:
        ax.set_ylabel(iso_y, color=positions.loc[iso_y]["color"])
    if ax_y == N_iso - 2:
        ax.set_xlabel(iso_x, color=positions.loc[iso_x]["color"])
plt.tight_layout()
plt.savefig(fig_fld / "dotplot.png", dpi=300)
plt.savefig(fig_fld / "dotplot.svg", dpi=300)
plt.show()
# %%

from Bio import Phylo
import seaborn as sns

tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, format="newick")

subtr = tree.root[1][1][1][1][1][1]
# subtr = tree.root[1][1][1][1][1][1][3]
isos = [x.name for x in subtr.get_terminals()]

fig, ax = plt.subplots(1, 1, figsize=(4, 6))
Phylo.draw(
    subtr,
    axes=ax,
    do_show=False,
    label_func=lambda x: "",
)

for n, leaf in enumerate(subtr.get_terminals()):
    iso = leaf.name
    if iso in positions.index:
        x = subtr.distance(leaf) + subtr.branch_length
        y = n + 1
        color = positions.loc[iso]["color"]
        plt.plot(x, y, marker="o", color=color)
        plt.text(3.2e-5, y, iso, color=color, ha="left", va="center", fontsize=12)

plt.title("")
plt.ylabel("")
plt.yticks([])
plt.xlabel("")
sns.despine(left=True)
plt.tight_layout()
plt.savefig(fig_fld / "tree.png")
plt.savefig(fig_fld / "tree.svg")
plt.show()

# %%
