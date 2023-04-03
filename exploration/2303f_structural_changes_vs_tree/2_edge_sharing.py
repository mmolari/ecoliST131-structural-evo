# %%

from utils import *

import pathlib
import copy

import pypangraph as pp
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import Phylo
from collections import defaultdict, Counter
from functools import cache

fig_fld = pathlib.Path("figs")
fig_fld.mkdir(exist_ok=True)


def svfig(svname):
    for suffix in ["png", "pdf"]:
        plt.savefig(fig_fld / f"{svname}.{suffix}", dpi=300)


# %%

# load graph and tree
pangraph_file = "../../results/ST131/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)

df = pan.to_blockstats_df()
core = df["core"].to_dict()
Ls = df["len"].to_dict()

tree = Phylo.read(
    "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk", format="newick"
)
tree.ladderize()

# extract path and edge representation
strains = set(pan.strains())
N = len(strains)
paths = {k: path_to_blocklist(pan, k) for k in strains}
edges = {k: blocklist_to_edgelist(bl) for k, bl in paths.items()}
# %%
#  extract edges
edge_strains = defaultdict(set)
for k, el in edges.items():
    for e in el:
        edge_strains[e].add(k)
print("n. unique edges:", len(edge_strains))
print("n. unique blocks:", len(pan.blocks))
# %%

# histogram of sharing of edges (might be similar to blocks, but more information)
L = [len(s) for k, s in edge_strains.items()]

kwargs = {
    "histtype": "step",
    "cumulative": True,
    "density": True,
    "bins": np.arange(N + 2) - 0.5,
}
plt.hist(L, label="edges", **kwargs)
plt.hist(df["n. strains"], label="blocks", **kwargs)
plt.legend(loc="upper left")
plt.xlabel("n. isolates")
plt.ylabel("cumulative distr.")
sns.despine()
plt.tight_layout()
svfig("edge_frequency")
plt.show()

# %%

# select edges that are present in only two strains
only_2 = [tuple(s) for e, s in edge_strains.items() if len(s) == 2]
missing_2 = [tuple(strains - s) for e, s in edge_strains.items() if len(s) == N - 2]
on = Counter(only_2).most_common()
ms = Counter(missing_2).most_common()
on = Counter(missing_2 + only_2).most_common()
# %%

tree_pairs = [n for n in tree.get_nonterminals() if len(n.get_terminals()) == 2]
tree_pairs = [tuple(set(x.name for x in n.get_terminals())) for n in tree_pairs]

tdfo = pd.DataFrame(
    {
        "p1": [p[0] for p, n in on],
        "p2": [p[1] for p, n in on],
        "n": [n for p, n in on],
        "on_tree": [p in tree_pairs for p, n in on],
    }
)
tdfo["ppv"] = np.cumsum(tdfo["on_tree"]) / np.cumsum(np.ones_like(tdfo["on_tree"]))

# %%
block_splits = pangraph_to_block_splits(pan)
block_pairs = [tuple(st) for b, st in block_splits.items() if len(st) == 2]
block_pairs = Counter(block_pairs).most_common()

tdfb = pd.DataFrame(
    {
        "p1": [p[0] for p, n in block_pairs],
        "p2": [p[1] for p, n in block_pairs],
        "n": [n for p, n in block_pairs],
        "on_tree": [p in tree_pairs for p, n in block_pairs],
    }
)
tdfb["ppv"] = np.cumsum(tdfb["on_tree"]) / np.cumsum(np.ones_like(tdfb["on_tree"]))


# %%

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(8, 6))

rank = np.arange(1, len(tdfo) + 1)

ax = axs[0]
ax.plot(rank, tdfo["n"], "o")
ax.set_ylabel("edge pair support")


ax = axs[1]
ax.plot(rank, tdfo["ppv"], "C0-", label="edge support")
ax.set_ylabel("p.p.v.")


rank = np.arange(1, len(tdfb) + 1)

ax = axs[1]
ax.plot(rank, tdfb["ppv"], "C1-", label="block support")
ax.axhline(len(tree_pairs) / (N * (N - 1) * 0.5), ls="--", label="chance", color="gray")
ax.set_ylim(bottom=0)
ax.legend()

ax = axs[2]
ax.plot(rank, tdfb["n"], "C1o")
ax.set_ylabel("block pair support")

ax.set_xscale("log")
ax.set_xlabel("rank of pair of isolates")
sns.despine()
plt.tight_layout()
svfig("edge_vs_block_support")
plt.show()
# %%


@cache
def mrca_clade_size(tree, ls):
    ns = [n for n in tree.get_terminals() if n.name in ls]
    mrca = tree.common_ancestor(ns)
    return len(mrca.get_terminals())


from itertools import combinations

all_mrcas = [mrca_clade_size(tree, ls) for ls in combinations(list(strains), 2)]
# %%

tdfb["mrca"] = tdfb.apply(lambda r: mrca_clade_size(tree, (r["p1"], r["p2"])), axis=1)
tdfo["mrca"] = tdfo.apply(lambda r: mrca_clade_size(tree, (r["p1"], r["p2"])), axis=1)

# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 4), sharex=True)
bins = np.arange(1, N + 2)

top = 20

kwargs = {"density": True, "bins": bins, "histtype": "step", "cumulative": True}
ax.hist(tdfo["mrca"].iloc[:top], **kwargs, label=f"top {top} edge ranking")
ax.hist(tdfb["mrca"].iloc[:top], **kwargs, label=f"top {top} block ranking")
ax.hist(all_mrcas, **kwargs, color="gray", ls="--", label="chance")
ax.set_xscale("log")
ax.set_xlim(left=1)
ax.set_xlabel("mrca clade size")
ax.set_ylabel("cumulative distr.")
ax.legend(loc="upper left")
sns.despine()
plt.tight_layout()
svfig(f"mrca_clade_size_top_{top}")
plt.show()

# %%

isolated_clade = {
    "NZ_CP076689",
    "NZ_CP076687",
    "NZ_CP014495",
    "NZ_JAOSEP010000001",
    "NZ_JAOSEH010000001",
    "NZ_CP014522",
    # "NZ_JAOSFC010000001"
}
complementary_clade = strains - isolated_clade

bc_e_support = [
    e
    for e, s in edge_strains.items()
    if (s == isolated_clade) or (s == complementary_clade)
]
bc_b_support = [
    b
    for b, s in block_splits.items()
    if (s == isolated_clade) or (s == complementary_clade)
]
print(len(bc_e_support))
print(len(bc_b_support))


# %%
