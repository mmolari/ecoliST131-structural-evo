# %%

# looking for single accessory blocks always flanked by core blocks

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import utils as ut

from collections import defaultdict
from Bio import Phylo

svfld = ut.fig_fld / "simple_flanking_dupl"
svfld.mkdir(exist_ok=True, parents=True)

# %%

pan = ut.load_pangraph()
paths = ut.pangraph_to_path_dict(pan)
bdf = pan.to_blockstats_df()
is_core = bdf["core"]

# %%
Js = {iso: ut.to_junctions(path.nodes) for iso, path in paths.items()}

# %%

# look for non-duplicated accessory blocks
mask = (~bdf["core"]) & bdf["duplicated"]
candidates = bdf[mask].index.to_numpy()
print(f"n. accessory candidates = {len(candidates)}")

# %%

# check if they are flanked by core blocks
accepted = list(candidates)

for iso, J in Js.items():
    for j in J:
        bid = j.center.id
        if not bid in accepted:
            continue
        l, r = j.left.id, j.right.id
        if not is_core[l] or not is_core[r]:
            accepted.remove(bid)

print(len(accepted))
# %%
