# %%

import numpy as np
import pandas as pd
import pathlib
from Bio import Phylo
import pypangraph as pp
import matplotlib.pyplot as plt

import utils as ut

from collections import defaultdict, Counter

svfld = pathlib.Path("figs/2b")
svfld.mkdir(parents=True, exist_ok=True)

# %%

prefix = "../../results/ST131/pangraph"
tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

# load tree and pangraph
tree = Phylo.read(tree_file, "newick")
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
is_core = bdf["core"].to_dict()
B_core = bdf[bdf["core"]].index.to_list()
B_len = bdf[bdf["core"]]["len"].to_dict()
N = len(pan.paths)

# %%
B_std = pd.read_csv("data/block_order.csv")

# %%

# 'NZ_CP104846': 2,
# 'NZ_CP104848': 2,
# 'NZ_JAOSCA010000001': 2,
# 'NZ_CP049085': 2,
# 'NZ_JAOSBZ010000001': 3,
# 'NZ_SEVU01000007': 2,
# 'NZ_SEVJ01000001': 4,
# 'NZ_JAOSEU010000001': 2,
# 'NZ_CP035377': 6,
# 'NZ_SEVM01000001': 2,
# 'NZ_JAOSEC010000001': 2,
# 'NZ_JAOSCQ010000001': 2,
# 'NZ_JAOSCB010000001': 2,


def dotplot(iso1, iso2, ax, roll=0, l_const=False):
    p1 = pan.paths[iso1]
    p2 = pan.paths[iso2]

    def core_path(p):
        B = p.block_ids
        S = p.block_strands
        B = np.roll(B, roll)
        S = np.roll(S, roll)
        return [(b, s) for b, s in zip(B, S) if is_core[b]]

    C1 = core_path(p1)
    C2 = core_path(p2)

    def coord_dict(C):
        x = 0
        coord = {}
        for b, s in C:
            l = 1 if l_const else B_len[b]
            coord[b] = (x, x + l) if s else (x + l, x)
            x += l
        return coord

    X1 = coord_dict(C1)
    X2 = coord_dict(C2)

    for b, s in C1:
        x = X1[b]
        y = X2[b]

        ax.plot(x, y)

    ax.set_xlabel(iso1)
    ax.set_ylabel(iso2)


iso1 = "NZ_CP035477"
iso2 = "NZ_JAOSCQ010000001"

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
dotplot(iso1, iso2, ax, l_const=False, roll=0)
# ax.set_xlim(1.3e6, 1.5e6)
# ax.set_ylim(0e6, 0.2e6)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(svfld / f"dotplot_{iso1}__{iso2}.png")
plt.show()
# %%

# fraction of terminal branches vs internal branches
L = tree.total_branch_length()
I = sum([n.branch_length for n in tree.get_terminals()])

# p-value for many rare edges on terminal branches
import scipy.stats as sps

sps.binom_test(11, 12, I / L)

# %%
