# %%
import json

import pypangraph as pp
import pandas as pd

from collections import defaultdict

import utils as ut

# %%

# load pangraph
pan = pp.Pangraph.load_json("data/pangraph.json")
bdf = pan.to_blockstats_df()
is_core = bdf["core"].to_dict()
block_len = bdf["len"].to_dict()
strains = pan.strains()
N_iso = len(strains)

# %%

# backbone blocks and edges
df = pd.read_csv("data/core_edge_frequencies.csv")
df = df[df["count"] == N_iso]
backbone_edges = [ut.Edge.from_str_id(e) for e in df["edge"]]
backbone_blocks = set(
    [e.left.id for e in backbone_edges] + [e.right.id for e in backbone_edges]
)

# %%

edges_pos = {}
for e in backbone_edges:
    if e.left.id > e.right.id:
        e = e.invert()
    assert e.left.id < e.right.id
    ln, rn = e.left, e.right
    pos_l, pos_r = pan.blocks[ln.id].alignment.pos, pan.blocks[rn.id].alignment.pos
    e_key = e.to_str_id()
    edges_pos[e_key] = defaultdict(lambda: [None, None, None])

    for aln_key in pos_l:
        iso, occ, strand = aln_key
        b, e = pos_l[aln_key]
        assert occ == 1
        same_strand = ln.strand == strand
        start = b if same_strand else e
        edges_pos[e_key][iso][0] = start
        edges_pos[e_key][iso][2] = same_strand

    for aln_key in pos_r:
        iso, occ, strand = aln_key
        b, e = pos_r[aln_key]
        assert occ == 1
        same_strand = rn.strand == strand
        end = e if same_strand else b
        edges_pos[e_key][iso][1] = end
        assert edges_pos[e_key][iso][2] == same_strand

    edges_pos[e_key] = dict(edges_pos[e_key])


# %%
