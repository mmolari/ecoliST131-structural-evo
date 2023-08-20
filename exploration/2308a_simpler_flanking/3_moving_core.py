# %%

import numpy as np
import pandas as pd

import utils as ut

from collections import Counter, defaultdict

LEN_THR = 500

# %%
pan = ut.load_pangraph()
paths = ut.pangraph_to_path_dict(pan)
bdf = pan.to_blockstats_df()
is_core = bdf["core"].to_dict()
block_len = bdf["len"].to_dict()


def core_filter(bid):
    if not is_core[bid]:
        return False
    if block_len[bid] < LEN_THR:
        return False
    return True


core_paths = ut.filter_paths(paths, core_filter)
strains = pan.strains()
N_iso = len(strains)
core_blocks = set(bdf[bdf["core"]].index.to_numpy())
# %%

# core-edges counter
Edges = []
for iso, path in core_paths.items():
    next_path = np.roll(path.nodes, -1)
    for l, r in zip(path.nodes, next_path):
        edge = ut.Edge(l, r)
        Edges.append(edge)
Edges = Counter(Edges)

backbone_edges = set([e for e, n in Edges.items() if n == N_iso])

with open(ut.expl_fld / "backbone_edges.txt", "w") as f:
    for e in backbone_edges:
        f.write(f"{e.to_str_id()}\n")
# %%

backbone_paths = []
path_idx_dict = {}
for e in backbone_edges:
    # print(f"-----\nedge = {e}")
    lid, rid = e.left.id, e.right.id
    left_assigned = lid in path_idx_dict
    right_assigned = rid in path_idx_dict

    if left_assigned and right_assigned:
        l_pid = path_idx_dict[lid]
        r_pid = path_idx_dict[rid]
        # if both belong to same path, circularization
        if l_pid == r_pid:
            # print("circularization")
            continue
        # else, join paths
        l_path = backbone_paths[l_pid]
        r_path = backbone_paths[r_pid]

        # print(f"join paths {l_path} and {r_path} with e={e}")

        if l_path.nodes[-1].id != lid:
            l_path = l_path.invert()
            assert l_path.nodes[-1].id == lid, f"path = {l_path}, e = {e}"
        if r_path.nodes[0].id != rid:
            r_path = r_path.invert()
            assert r_path.nodes[0].id == rid, f"path = {r_path}, e = {e}"
        assert l_path.nodes[-1].strand == e.left.strand, f"path = {l_path}, e = {e}"
        assert r_path.nodes[0].strand == e.right.strand, f"path = {r_path}, e = {e}"

        new_path = ut.Path(l_path.nodes + r_path.nodes)
        backbone_paths[l_pid] = None
        backbone_paths[r_pid] = None
        backbone_paths.append(new_path)
        # print(f"added path n. {len(backbone_paths) - 1}")
        for n in new_path.nodes:
            path_idx_dict[n.id] = len(backbone_paths) - 1

        # print(f"outcome: {new_path}")

    elif right_assigned:
        # append to left of path
        pid = path_idx_dict[rid]
        path = backbone_paths[pid]
        # print(f"right-assigned, e={e}, p={path}")

        if path.nodes[0].id != rid:
            path = path.invert()
            assert path.nodes[0].id == rid, f"path = {path}, e = {e}"
        assert path.nodes[0].strand == e.right.strand, f"path = {path}, e = {e}"
        path.nodes.insert(0, e.left)
        path_idx_dict[lid] = pid
        backbone_paths[pid] = path
        # print(f"outcome: {path}")
    elif left_assigned:
        # append to right of path
        pid = path_idx_dict[lid]
        path = backbone_paths[pid]
        # print(f"left-assigned, e={e}, p={path}")

        if path.nodes[-1].id != lid:
            path = path.invert()
            assert path.nodes[-1].id == lid, f"path = {path}, e = {e}"
        assert path.nodes[-1].strand == e.left.strand, f"path = {path}, e = {e}"
        path.nodes.append(e.right)
        path_idx_dict[rid] = pid
        backbone_paths[pid] = path
        # print(f"outcome: {path}")
    else:
        # create new path
        path = ut.Path([e.left, e.right])
        backbone_paths.append(path)
        path_idx_dict[lid] = len(backbone_paths) - 1
        path_idx_dict[rid] = len(backbone_paths) - 1
        # print(f"new path n. {len(backbone_paths) - 1}")

backbone_paths = [p for p in backbone_paths if p is not None]
backbone_paths
# %%
