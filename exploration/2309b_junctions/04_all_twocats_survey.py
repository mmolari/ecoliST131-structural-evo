# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import numpy as np
import os
import json

import pypangraph as pp
import utils as ut
import mugration_utils as mu


from collections import defaultdict
from Bio import Phylo


fig_fld = pathlib.Path("figs")


fld = pathlib.Path("../../results/ST131")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

# %%
df = pd.read_csv(df_file, index_col=0)
Js = df.index.to_list()
# %%

tree = Phylo.read(tree_file, "newick")
tree.ladderize()
branch_len = {b.name: b.branch_length for b in tree.get_terminals()}
branch_len |= {b.name: b.branch_length for b in tree.get_nonterminals()}

# %%


def event_category(p1, p2):
    p1, p2 = set(p1), set(p2)
    if p1.issubset(p2):
        return "gain", list(p2 - p1)
    elif p2.issubset(p1):
        return "loss", list(p1 - p2)
    else:
        return "other", list((p1 | p2) - (p1 & p2))


def are_clade(tree, isolates):
    "check if the isolates form a clade in the tree"
    # get common ancestor
    leaves = {leaf.name: leaf for leaf in tree.get_terminals()}
    common_ancestor = tree.common_ancestor(*[leaves[i] for i in isolates])
    # check if terminal nodes are equal to clade
    clade = set([leaf.name for leaf in common_ancestor.get_terminals()])
    if clade == set(isolates):
        return common_ancestor.name
    else:
        return None


# def find_coordinates(pan, iso, block):
# path = pan.paths[iso]


events_df = {}
for j in Js:
    info = {}
    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    paths = ut.pangraph_to_path_dict(pan)
    path_cat = ut.path_categories(paths)

    if not len(path_cat) == 2:
        continue

    n1, p1, i1 = path_cat[0]
    n2, p2, i2 = path_cat[1]

    iso = None
    if n2 == 1:
        info["type"] = "terminal"
        iso = i2[0]
    else:
        info["type"] = "internal"
        iso = are_clade(tree, i2)

    if not (iso is None):
        info["isolate"] = iso
        info["branch len"] = branch_len[iso]

    cat, blocks = event_category(p1, p2)
    info["event category"] = cat
    if len(blocks) == 1:
        bid = blocks[0].id
        info["block"] = bid
        bdf = pan.to_blockstats_df()
        info["block len"] = bdf.loc[bid, "len"]

    events_df[j] = info

events_df = pd.DataFrame.from_dict(events_df, orient="index")
events_df.to_csv("data/twocat_events.csv")
events_df
# %%

sns.histplot(events_df, x="event category", hue="type", multiple="stack")
plt.show()

# %%
sns.histplot(
    events_df,
    x="block len",
    hue="event category",
    log_scale=True,
    stat="probability",
    common_norm=False,
    common_bins=True,
    # multiple="stack",
)
plt.show()

# %%
