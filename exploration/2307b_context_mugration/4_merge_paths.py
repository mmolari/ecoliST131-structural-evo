# %%
import numpy as np
import pandas as pd
import utils as ut

from collections import Counter

pan = ut.load_pangraph()
bdf = pan.to_blockstats_df()
# %%

# only keep selected blocks
retained_blocks = np.loadtxt(ut.expl_fld / "retained_blocks.txt", dtype=str)
bdf = bdf.loc[retained_blocks].copy()
block_PA = pan.to_blockcount_df().T.loc[bdf.index] > 0

paths = ut.pangraph_to_path_dict(pan)
paths = ut.filter_paths(paths, lambda bid: bid in retained_blocks)

# %%


def path_to_edges(path):
    edges = []
    L = len(path.nodes)
    for l in range(L):
        left = path.nodes[(l - 1) % L]
        right = path.nodes[l]
        edge = ut.Edge(left, right)
        edges.append(edge)
    return edges


iso_edges = {iso: path_to_edges(path) for iso, path in paths.items()}
edge_count = Counter(sum(iso_edges.values(), start=[]))

merge_list = {}
for e in edge_count.keys():
    edge_pa = {iso: (e in E) for iso, E in iso_edges.items()}
    left, right = e.left.id, e.right.id
    # same PA pattern of edge and block
    same_left = np.all(block_PA.loc[left] == pd.Series(edge_pa))
    same_right = np.all(block_PA.loc[right] == pd.Series(edge_pa))

    if not (same_left and same_right):
        continue

    if (left in merge_list) and (right in merge_list):
        ref = merge_list[left]
        old_ref = merge_list[right]
        for k in merge_list:
            if merge_list[k] == old_ref:
                merge_list[k] = ref
    elif left in merge_list:
        ref = merge_list[left]
    elif right in merge_list:
        ref = merge_list[right]
    else:
        ref = left
    merge_list[left] = ref
    merge_list[right] = ref


# %%

# clean up paths
kept = set(merge_list.values())
deleted = set(merge_list.keys()) - kept
for iso, path in paths.items():
    for l in range(len(path.nodes))[::-1]:
        if path.nodes[l].id in deleted:
            path.nodes.pop(l)

# adapt stats
bdf["parts"] = 1
for k in deleted:
    v = merge_list[k]
    bdf.loc[v, "len"] += bdf.loc[k, "len"]
    bdf.loc[v, "parts"] += 1
    assert bdf.loc[v, "count"] == bdf.loc[k, "count"]
    bdf.drop(k, inplace=True)

# %%
