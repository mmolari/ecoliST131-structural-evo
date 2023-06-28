# %%

import numpy as np
import pandas as pd
import pathlib
from Bio import Phylo
import pypangraph as pp

import utils as ut

from collections import defaultdict, Counter


svfld = pathlib.Path("data")
svfld.mkdir(exist_ok=True, parents=True)

# initialize dataframe
df = {}
# %%

prefix = "../../results/ST131/pangraph"
tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

# load tree and pangraph
tree = Phylo.read(tree_file, "newick")
pan = pp.Pangraph.load_json(pangraph_file)

# %%

# generate block presence/absence matrix
df = pan.to_blockcount_df()

# select only non-duplicated genes
dupl_mask = (df > 1).any(axis=0)
df = df.loc[:, ~dupl_mask]
# %%

# infer ancestral block presence-absence state
B = df.columns.to_list()

PA = {}

for clade in tree.find_clades(order="postorder"):
    if clade.is_terminal():
        PA[clade.name] = (df.loc[clade.name, :] > 0).to_dict()
    else:
        C = clade.clades
        sdf = pd.DataFrame({c.name: PA[c.name] for c in C}).T
        # only infer as present if present in all children
        PA[clade.name] = sdf.all(axis=0).to_dict()


# %%

# ideas: what if one updates the tree with maximum parsimony?
