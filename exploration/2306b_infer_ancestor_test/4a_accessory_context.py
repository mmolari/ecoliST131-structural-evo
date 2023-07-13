# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
import utils4 as ut

from Bio import Phylo
from collections import defaultdict, Counter


# %%

# load pangraph
pan = ut.load_pangraph()
bdf = pan.to_blockstats_df()
is_core = bdf["core"].to_dict()
is_dupl = bdf["duplicated"].to_dict()
bl_Ls = bdf["len"].to_dict()
# %%

LEN_THRESH = 500


def filter_func(bid):
    if bl_Ls[bid] < LEN_THRESH:
        return False
    if is_dupl[bid]:
        return False
    return True


def to_nodelist(path):
    B = path.block_ids
    S = path.block_strands
    nodes = [ut.Node(b, s) for b, s in zip(B, S)]
    return nodes


def to_acc_adjacencies(nodes):
    adj = {}
    N = len(nodes)

    # appen initial blocks to the end of the list for periodic boundary conditions
    nodes = [nodes[-1]] + list(nodes) + [nodes[0]]

    for i in range(1, N + 1):
        n = nodes[i]
        if is_core[n.id]:
            continue
        adj[n.id] = ut.Junction(nodes[i - 1], n, nodes[i + 1])

    return adj


def to_core_adjacencies(nodes):
    adj = {}
    N = len(nodes)

    for i, n in enumerate(nodes):
        if is_core[n.id]:
            continue
        bef, aft = None, None
        bi, ai = i, i
        while bef is None:
            bi = (bi - 1) % N
            if is_core[nodes[bi].id]:
                bef = nodes[bi]
        while aft is None:
            ai = (ai + 1) % N
            if is_core[nodes[ai].id]:
                aft = nodes[ai]
        adj[n.id] = ut.Junction(bef, n, aft)

    return adj


Adj = {}

for path in pan.paths:
    name = path.name
    print(f"processing {name}")
    nodes = to_nodelist(path)
    print(f"  {len(nodes)} blocks")

    # filter out short and duplicated blocks
    mask = np.array([filter_func(n.id) for n in nodes])
    nodes = np.array(nodes)[mask]
    print(f"  {np.sum(~mask)} blocks filtered")
    print(f"  {len(nodes)} blocks remaining")

    adj_acc = to_acc_adjacencies(nodes)
    adj_core = to_core_adjacencies(nodes)

    Adj[name] = {"acc": adj_acc, "core": adj_core}

    # filter paths

# %%

ctr_a, ctr_c = defaultdict(Counter), defaultdict(Counter)
for iso in pan.strains():
    adj_core, adj_acc = Adj[iso]["core"], Adj[iso]["acc"]
    for k, j in adj_core.items():
        ctr_c[k].update([j])
    for k, j in adj_acc.items():
        ctr_a[k].update([j])

ctr_df = []
for k in ctr_c:
    c_ct = ctr_c[k]
    a_ct = ctr_a[k]
    ctr_df.append(
        {
            "block": k,
            "core_ctx": len(c_ct),
            "acc_ctx": len(a_ct),
            "N": sum(c_ct.values()),
        }
    )
ctr_df = pd.DataFrame(ctr_df)
ctr_df

# %%

sns.histplot(data=ctr_df, x="core_ctx", y="acc_ctx", bins=[np.arange(12) + 0.5] * 2)

# %%
mask = ctr_df["core_ctx"] == ctr_df["acc_ctx"]
ctr_df[mask]

# %%
