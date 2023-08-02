# %%

import json

import pandas as pd
import numpy as np
import utils as ut

from collections import defaultdict

sv_fld = ut.expl_fld / "filtered_paths"

pan = ut.load_pangraph()
bdf = pd.read_csv(sv_fld / "block_stats.csv", index_col=0)
is_core = bdf["core"].to_dict()

# exclude core
non_core_bids = bdf[~bdf["core"]].index.to_numpy()

PA = pan.to_blockcount_df()
PA = PA[non_core_bids] > 0

# %% 3. extract block PA
for bid, col in PA.iteritems():
    PA[bid] = PA[bid].map({True: "P", False: "-"})
PA.index.name = "#name"
PA.to_csv(sv_fld / "PA_simple.csv")

# %%
# load paths
with open(sv_fld / "paths.json", "r") as f:
    paths = json.load(f)
paths = {iso: ut.Path.from_list(p) for iso, p in paths.items()}
# %%

path_adj = {}
for iso, p in paths.items():
    path_adj[iso] = ut.to_core_adjacencies(p.nodes, is_core)
path_adj


PA_ctx = PA.copy()
PA_ctx.loc[:, :] = "-"

adj_letter = defaultdict(lambda: {"next": "A"})
for iso, path in path_adj.items():
    for j in path:
        b = j.center.id
        if not (j in adj_letter[b]):
            adj_letter[b][j] = adj_letter[b]["next"]
            next_letter = chr(ord(adj_letter[b]["next"]) + 1)
            adj_letter[b]["next"] = next_letter
        PA_ctx.loc[iso, b] = adj_letter[b][j]

PA_ctx.index.name = "#name"
PA_ctx.to_csv(sv_fld / "PA_context.csv")

PA_ctx
# %%
