# %%

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
from Bio import Phylo
import pypangraph as pp
import seaborn as sns

import utils as ut

from collections import defaultdict, Counter


# %%

svfld = pathlib.Path("figs/1a")
svfld.mkdir(exist_ok=True, parents=True)


def svfig(name):
    plt.savefig(svfld / name, dpi=300, facecolor="white", bbox_inches="tight")


# %%

prefix = "../../results/ST131/pangraph"
tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

# load tree
tree = Phylo.read(tree_file, "newick")
# load pangraph
pan = pp.Pangraph.load_json(pangraph_file)

# %%

# get terminal branch lengths
len_tips = {}
for t in tree.get_terminals():
    len_tips[t.name] = t.branch_length
isolates = list(len_tips.keys())

# %%

bpa = pan.to_blockcount_df()
bdf = pan.to_blockstats_df()

# select only terminal blocks
mask = (bpa > 0).sum() == 1
t_blocks = mask.index[mask].tolist()

n_acc = defaultdict(int)  # number of acc blocks
n_dupl = defaultdict(int)  # number of duplicated blocks
L_acc = defaultdict(int)  # total length of acc blocks
L_dupl = defaultdict(int)  # total length of duplicated blocks
acc_list = defaultdict(list)  # list of acc blocks per isolate (non-duplicated)
for i in isolates:
    for b in t_blocks:
        if bpa.loc[i, b] > 0:
            if bdf["duplicated"][b]:
                n_dupl[i] += 1
                L_dupl[i] += bdf["len"][b]
            else:
                acc_list[i].append(b)
                n_acc[i] += 1
                L_acc[i] += bdf["len"][b]


def block_groups(path, selected_blocks):
    """path is a list of blocks, with circular boundary conditions.
    This function groups adjacent selected blocks into groups.
    """
    p = np.array(path)

    s = np.isin(p, selected_blocks).astype(int)
    breakpoints = np.abs(s - np.roll(s, 1))
    return np.sum(breakpoints) / 2


n_events = {}
for i, B in acc_list.items():
    p = pan.paths[i]
    bl_list = [b for b in p.block_ids if not bdf["duplicated"][b]]
    n_events[i] = block_groups(bl_list, B)

# %%

df = pd.DataFrame(
    {
        "branch_len": len_tips,
        "n_acc": n_acc,
        "n_dupl": n_dupl,
        "L_acc": L_acc,
        "L_dupl": L_dupl,
        "n_events": n_events,
    }
)
df
# %%
sns.histplot(data=df, x="branch_len", y="n_acc", bins=20)
plt.show()

sns.histplot(data=df, x="branch_len", y="L_acc", bins=20)
plt.show()

sns.histplot(data=df, x="branch_len", y="n_dupl", bins=20)
plt.show()

sns.histplot(data=df, x="branch_len", y="L_dupl", bins=20)
plt.show()

sns.histplot(data=df, x="branch_len", y="n_events", bins=20)
plt.show()

sns.histplot(data=df, x="n_acc", y="n_events", bins=20)
plt.show()
# %%

# duplicated genes: instead of presence/absence look at junctions

junction_list = defaultdict(list)
for i in isolates:
    p = pan.paths[i]
    blks = p.block_ids
    blks_left = np.roll(blks, 1)
    blks_right = np.roll(blks, -1)
    strands = p.block_strands
    strands_left = np.roll(strands, 1)
    strands_right = np.roll(strands, -1)

    for k in range(len(blks)):
        if bdf["duplicated"][blks[k]]:
            n_l = ut.Node(blks_left[k], strands_left[k])
            n_r = ut.Node(blks_right[k], strands_right[k])
            n_c = ut.Node(blks[k], strands[k])
            j = ut.Junction(n_l, n_r, n_c)
            junction_list[i].append(j)

# %%

# find junctions present in only one isolate
j_counts = Counter([j for i in junction_list for j in junction_list[i]])

# TODO
# relate this to terminal branches
# then write this into a story

# %%
