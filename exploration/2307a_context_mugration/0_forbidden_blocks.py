# %%
"""
Find forbidden blocks:
accessory blocks that are found at least once in a core-genome synteny breakpoint.
"""

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


# %% 2a: clean paths

THR_LEN = 500

pan = ut.load_pangraph()
bdf = pan.to_blockstats_df()
block_Ls = bdf["len"].to_dict()
is_dupl = bdf["duplicated"].to_dict()
is_core = bdf["core"].to_dict()

# core-genes paths
paths = {}
for p in pan.paths:
    B = p.block_ids
    S = p.block_strands
    path = []
    for b, s in zip(B, S):
        if block_Ls[b] < THR_LEN:
            continue
        if not is_core[b]:
            continue
        node = ut.Node(b, s)
        path.append(node)
    paths[p.name] = path

# %%
Edges = []

for iso, path in paths.items():
    next_path = np.roll(path, -1)
    for l, r in zip(path, next_path):
        edge = ut.Edge(l, r)
        Edges.append(edge)
Edges = Counter(Edges)

# %%
N = len(ut.load_pangraph().strains())
forbidden_edges = set([e for e, n in Edges.items() if n < N])

# %%

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

# %%
path_adj = {}
for iso, p in paths.items():
    path_adj[iso] = ut.to_core_adjacencies(p, is_core)

# %%

forbidden_blocks = []
for iso, ap in path_adj.items():
    for bid, j in ap.items():
        if is_core[bid]:
            continue
        edge = ut.Edge(j.left, j.right)
        if edge in forbidden_edges:
            forbidden_blocks.append(bid)
forbidden_blocks = set(forbidden_blocks)

# %%

# save list of forbidden blocks
with open(ut.forb_blocks_file, "w") as f:
    for bid in forbidden_blocks:
        f.write(f"{bid}\n")

# %%
