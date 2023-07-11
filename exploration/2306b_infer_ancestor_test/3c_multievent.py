# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import json

import pypangraph as pp
import utils3 as ut

from Bio import Phylo
from collections import Counter

fld = pathlib.Path("data/mugration")
fld.mkdir(exist_ok=True)

fig_fld = pathlib.Path("figs/3c")
fig_fld.mkdir(exist_ok=True)

# %%

prefix = "../../results/ST131/pangraph"
original_tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
tree_file = fld / "tree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

# load tree and pangraph
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()

# load gene presence/absence inference
pa_inference_file = fld / "infer_pa.json"
with open(pa_inference_file, "r") as f:
    pa_inference = json.load(f)


# %%

# look at patterns for multiple events for long blocks
df = pd.read_csv(fld / "infer_pa_summary.csv", index_col=0)

mask = df["category"] == "multiple-events"
sdf = df[mask].sort_values("len", ascending=False)
sdf

# %%
mask = sdf["len"] > 10000
for bid in sdf[mask].index.to_list():
    fig, ax = ut.plot_tree_events(tree_file, pa_inference, bid)
    ax.set_title(f"block {bid} - len {sdf.loc[bid, 'len']//1000} kbp")
    plt.tight_layout()
    plt.savefig(fig_fld / f"tree_{bid}.png")
    plt.show()

# %%

# where was block JTBIYNWKTP integrated?
bid = "ZASKRVDPQW"
for p in pan.paths:
    B = p.block_ids
    S = p.block_strands
    if bid in B:
        cp = None
        ca = None
        met = False
        k = 0
        for b, s in zip(B, S):
            core = bdf.loc[b]["core"]
            if core:
                if not met:
                    cp = (b, s)
                else:
                    ca = (b, s)
                    break
            elif b == bid:
                print(b, s)
                print(B[k - 2 : k + 3], S[k - 2 : k + 3])
                met = True
            k += 1
        print(f"path {p.name} - core previous {cp} - core after {ca}")

# %%
