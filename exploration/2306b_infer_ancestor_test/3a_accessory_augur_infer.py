# %%

import pathlib
import os

import itertools as itt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json

import pypangraph as pp
import utils3 as ut

from Bio import Phylo
from collections import Counter

fld = pathlib.Path("data/mugration")
fld.mkdir(exist_ok=True)

# %%

prefix = "../../results/ST131/pangraph"
original_tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
tree_file = fld / "tree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

# load tree and pangraph
tree = Phylo.read(original_tree_file, "newick")
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()

for n, node in enumerate(tree.get_nonterminals()):
    node.name = f"node_{n:02d}"

Phylo.write(tree, tree_file, "newick", format_branch_length="%.5e")

# %%

# select non-duplicated accessory genes
mask = ~bdf["core"]
mask &= ~bdf["duplicated"]

# list of accessory blocks ids
B_acc = bdf[mask].index.to_list()
# %%

# extract presence/absence matrix
pa_file = fld / "pa_matrix.csv"
pa_df = pan.to_blockcount_df()
pa_df = pa_df.loc[:, B_acc]
pa_df.index.name = "#name"
pa_df.replace({0: "A", 1: "P"}).to_csv(pa_file)
# %%

pa_inference = {}

for n_attr, attr in enumerate(B_acc):
    print(f"processing {attr} \t - \t {n_attr+1}/{len(B_acc)}")
    out_dir = fld / attr

    # run treetime mugration
    cmd = f"""
    treetime mugration \
        --tree {tree_file} \
        --states {pa_file} \
        --attribute {attr} \
        --outdir {out_dir} \
    """

    os.system(cmd)

    # output file names
    nexus_file = out_dir / "annotated_tree.nexus"
    rates_file = out_dir / "GTR.txt"

    # parse rates
    rates = ut.parse_gtr(out_dir / "GTR.txt")

    # parse inferred presence/absence pattern
    pa_pattern = ut.parse_nexus(nexus_file)

    os.system(f"rm -r {out_dir}")

    pa_inference[attr] = {
        "pa_pattern": pa_pattern,
        "rates": rates,
    }

# write to json
pa_inference_file = fld / "infer_pa.json"
with open(pa_inference_file, "w") as f:
    json.dump(pa_inference, f)


# %%
