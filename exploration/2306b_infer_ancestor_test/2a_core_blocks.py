# %%

import numpy as np
import pandas as pd
import pathlib
from Bio import Phylo
import pypangraph as pp
import matplotlib.pyplot as plt

import utils as ut

from collections import defaultdict, Counter

svfld = pathlib.Path("figs/2a")
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
N = len(pan.paths)
# %%

# generate paths as list of nodes, only for core genes
paths = {}
for path in pan.paths:
    B = path.block_ids
    S = path.block_strands
    p = [ut.Node(B[i], S[i]) for i in range(len(path)) if is_core[B[i]]]
    paths[path.name] = p

# %%

junctions = {}
for iso, P in paths.items():
    L = len(P)
    J = {}
    for i in range(L):
        l = P[(i - 1) % L]
        c = P[i]
        r = P[(i + 1) % L]
        J[c.id] = ut.Junction(l, c, r)
    junctions[iso] = J

# %%

# find breakpoint blocks:
bkp = {}
for b in B_core:
    joints = defaultdict(list)
    for iso, J in junctions.items():
        joints[J[b]].append(iso)
    if len(joints) > 1:
        bkp[b] = joints
        print(f"block {b} -> {[len(v) for k,v in joints.items()]}")
bkp

# %%

edges = {}
for iso, P in paths.items():
    L = len(P)
    E = []
    for i in range(L):
        l = P[i]
        r = P[(i + 1) % L]
        E.append(ut.Edge(l, r))
    edges[iso] = E

# %%

iso_per_edge = defaultdict(list)
for iso, E in edges.items():
    for e in E:
        iso_per_edge[e].append(iso)

# %%

#  edge frequency distribution

edge_freq = {k: len(i) for k, i in iso_per_edge.items()}
plt.hist(edge_freq.values(), bins=np.arange(N + 2) - 0.5)
plt.xlabel("core-edge frequency")
plt.ylabel("n. edges")
plt.tight_layout()
plt.savefig(svfld / "core_edge_freq.png")
plt.show()
# %%
bkp_iso = []
for e, iso in iso_per_edge.items():
    if len(iso) == N:
        continue
    if len(iso) < N / 2:
        bkp_iso += iso
bkp_iso = Counter(bkp_iso)
bkp_iso
# 2 edges shared by 2 isolates
# 29 edges private to 1 isolate
# 13 isolates
# %%

fig, ax = plt.subplots(1, 1, figsize=(5, 10))
Phylo.draw(
    tree,
    label_func=lambda n: f"{n.name} - {bkp_iso[n.name]}" if n.name in bkp_iso else None,
    do_show=False,
    axes=ax,
)
plt.title("n. breakpoints in core genome order")
plt.tight_layout()
plt.savefig(svfld / "n_core_breakpoints.png")
plt.show()
# %%

# extract standard core block ordering
df = None
for iso in pan.strains():
    if iso not in bkp_iso:
        p = paths[iso]
        df = []
        for n in p:
            df.append({"block": n.id, "strand": n.strand})
        df = pd.DataFrame(df)
        break
df.to_csv("data/block_order.csv", index=False)

# %%
