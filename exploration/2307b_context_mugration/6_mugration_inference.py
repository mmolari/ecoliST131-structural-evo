# %%
import os
import pathlib
import json

import numpy as np
import pandas as pd

import utils as ut


loadfld = ut.expl_fld / "filtered_paths"

svfld = ut.expl_fld / "mugration"
svfld.mkdir(exist_ok=True)

# %%
bdf = pd.read_csv(loadfld / "block_stats.csv", index_col=0)
Bs = bdf[~bdf["core"]].index.to_numpy()

# %%


pa_inference = {}

for i, bid in enumerate(Bs):
    print(f"processing {bid} \t - \t {i+1}/{len(Bs)}")
    out_dir = svfld / bid

    # run treetime mugration
    states = loadfld / "PA_simple.csv"

    cmd = f"""
    treetime mugration \
        --tree {ut.named_nodes_tree_file} \
        --states {states} \
        --attribute {bid} \
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

    pa_inference[bid] = {
        "pa_pattern": pa_pattern,
        "rates": rates,
    }

# write to json
pa_inference_file = svfld / "infer_pa_simple.json"
with open(pa_inference_file, "w") as f:
    json.dump(pa_inference, f)

# %%


pa_inference = {}

for i, bid in enumerate(Bs):
    print(f"processing {bid} \t - \t {i+1}/{len(Bs)}")
    out_dir = svfld / bid

    # run treetime mugration
    states = loadfld / "PA_context.csv"

    cmd = f"""
    treetime mugration \
        --tree {ut.named_nodes_tree_file} \
        --states {states} \
        --attribute {bid} \
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

    pa_inference[bid] = {
        "pa_pattern": pa_pattern,
        "rates": rates,
    }

# write to json
pa_inference_file = svfld / "infer_pa_context.json"
with open(pa_inference_file, "w") as f:
    json.dump(pa_inference, f)

# %%
