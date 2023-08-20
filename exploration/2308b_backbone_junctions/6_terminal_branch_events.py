# %%

import pathlib

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pypangraph as pp
import utils as ut

from collections import defaultdict
from Bio import Phylo

fld = pathlib.Path("../../results/ST131/backbone_joints/asm20-100-5/")

# %%
df = pd.read_csv(
    fld / "junctions_stats.csv",
    index_col=0,
)
sdf = df[df["singleton"]]
Js = sdf.index.to_list()
del df, sdf

# %%

tree_file = "data/named_tree.nwk"
tree = Phylo.read(tree_file, "newick")
tree.ladderize()

terminal_len = {b.name: b.branch_length for b in tree.get_terminals()}

# %%


def path_categories(paths):
    """Returns a list of touples, one per non-empty path, with the following info:
    (count, path, [list of isolates])"""
    iso_list = defaultdict(list)
    n_paths = defaultdict(int)
    nodes = {}
    for iso, path in paths.items():
        if len(path.nodes) > 0:
            n_paths[path] += 1
            iso_list[path].append(iso)
            nodes[path] = path.nodes

    # sort by count
    path_cat = [(count, nodes[path], iso_list[path]) for path, count in n_paths.items()]
    path_cat.sort(key=lambda x: x[0], reverse=True)
    return path_cat


ev_ctr = defaultdict(int)
for j in Js:
    pan_file = fld / "joints_pangraph" / f"{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    paths = ut.pangraph_to_path_dict(pan)
    path_cat = path_categories(paths)

    assert len(path_cat) == 2
    assert path_cat[1][0] == 1
    iso = path_cat[1][2][0]
    ev_ctr[iso] += 1

# %%
df = pd.DataFrame.from_dict(terminal_len, orient="index", columns=["branch_length"])
df["ev_count"] = df.index.map(ev_ctr)
df["#events > 0"] = df["ev_count"] > 0
df
# %%
df["#events > 0"].value_counts()
# %%

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
ax = axs[0]
sns.scatterplot(data=df, x="branch_length", y="ev_count", alpha=0.5, ax=ax)
ax.set_xlabel("terminal branch length")
ax.set_ylabel("number of events")

ax = axs[1]
sns.histplot(data=df, x="branch_length", hue="#events > 0", element="step", ax=ax)
ax.set_xlabel("terminal branch length")
ax.set_ylabel("n. isolates")

sns.despine()
plt.tight_layout()
plt.savefig("figs/src_6/terminal_branch_events.png", facecolor="white")
plt.show()

# %%
