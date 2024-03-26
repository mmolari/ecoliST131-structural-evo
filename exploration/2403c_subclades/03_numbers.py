# %%
import pandas as pd
import numpy as np
from Bio import Phylo
import json

df = []

synt = {
    "ST131_ABC": 22,
    "ST131_sub_BC": 19,
    "ST131_sub_C": 19,
    "ST131_sub_C2": 14,
}

for dset in ["ST131_ABC", "ST131_sub_BC", "ST131_sub_C", "ST131_sub_C2"]:

    fname = f"../../results/{dset}/pangraph/asm20-100-5-filtered-coretree.nwk"
    tree = Phylo.read(fname, "newick")

    fname = f"../../results/{dset}/pangraph/asm20-100-5-alignment/filtered_corealignment_info_size.json"
    with open(fname) as f:
        aln_data = json.load(f)
    ca = aln_data["core aln size"]
    pa = aln_data["polished aln size"]

    fname = f"../../results/{dset}/rates/asm20-100-5/merged_events_df.csv"
    edf = pd.read_csv(fname, index_col=0)
    ntg = ((edf["type"] == "gain") & edf["terminal"]).sum()
    ntl = ((edf["type"] == "loss") & edf["terminal"]).sum()

    # tree length
    tree_len = tree.total_branch_length()
    tree_terminal_length = sum([t.branch_length for t in tree.get_terminals()])
    N = len(tree.get_terminals())
    df.append(
        {
            "dataset": dset,
            "tree_len": tree_len,
            "tree_terminal_len": tree_terminal_length,
            "n_iso": N,
            "noncons_synt_pattrn": synt[dset],
            "core aln size": ca,
            "polished aln size": pa,
            "synt. change rate (tree)": tree_len / synt[dset],
            "synt. change rate (red)": tree_len * pa / synt[dset],
            "gain ter": ntg,
            "loss ter": ntl,
            "gain ter rate": tree_terminal_length * pa / ntg,
            "loss ter rate": tree_terminal_length * pa / ntl,
        }
    )

df = pd.DataFrame(df).set_index("dataset")
df

# %%

# junction positions:

dset = "ST131_ABC"
# dset = "ST131_sub_BC"
# dset = "ST131_sub_C"
# dset = "ST131_sub_C2"

pos_file = f"../../results/{dset}/backbone_joints/asm20-100-5/joints_pos.json"
with open(pos_file) as f:
    jp = json.load(f)

# %%
iso = "NZ_CP107153.1"
p = 1739875
pos_iso = {k: v[iso] for k, v in jp.items() if iso in v}

for k, v in pos_iso.items():
    cb, ab, ae, ce, strand = v
    if cb <= p <= ce:
        print(k, v)

# %%
