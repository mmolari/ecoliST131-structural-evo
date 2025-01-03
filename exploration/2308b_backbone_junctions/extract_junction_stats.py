# %%
import pathlib

import numpy as np
import pandas as pd

import pypangraph as pp
import utils as ut

from scipy.stats import entropy
from collections import defaultdict, Counter

pan_fld = "../../results/ST131/backbone_joints/asm20-100-5/joints_pangraph/"

pan_files = list(pathlib.Path(pan_fld).glob("*.json"))
# %%

# entropy
# n.blocks
# n. path categories
# singletons
# min/max length
# left core length
# right core length


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


def get_stats(pan_file):
    # initialize stat dictionary
    stats = {"edge": pan_file.stem.split(".")[0]}

    pan = pp.Pangraph.load_json(pan_file)
    bdf = pan.to_blockstats_df()
    N = len(pan.strains())

    # number of blocks
    stats["n_blocks"] = len(bdf)

    # number of categories
    paths = ut.pangraph_to_path_dict(pan)
    path_cat = path_categories(paths)
    stats["n_categories"] = len(path_cat)
    stats["majority_category"] = path_cat[0][0]
    stats["singleton"] = path_cat[0][0] == N - 1
    stats["cat_entropy"] = entropy([count for count, _, _ in path_cat])

    # path lengths
    Ls = []
    len_entr = 0
    core_sides = {"left": [], "right": []}
    for count, nodes, isolates in path_cat:
        lengths = [bdf["len"][node.id] for node in nodes]
        Ls += [sum(lengths)] * count
        len_entr += entropy(lengths) * count
        core_sides["left"].append(nodes[0].id)
        core_sides["right"].append(nodes[-1].id)
    stats["min_length"] = min(Ls)
    stats["max_length"] = max(Ls)
    stats["mean_length"] = np.mean(Ls)
    stats["length_entropy"] = len_entr / N
    for side in ["left", "right"]:
        stats[f"core_{side}_length"] = None
        if not np.all(np.array(core_sides[side]) == core_sides[side][0]):
            continue
        if not bdf["core"][core_sides[side][0]]:
            continue
        stats[f"core_{side}_length"] = bdf["len"][core_sides[side][0]]

    return stats


df = []
for pan_file in pan_files:
    stats = get_stats(pan_file)
    df.append(stats)
df = pd.DataFrame(df)
df

# %%

# which non-empty blocks are in the backbone?
import seaborn as sns
import matplotlib.pyplot as plt

sns.histplot(data=df, x="n_blocks", bins=np.arange(100) - 0.5, hue="singleton")
plt.yscale("log")
plt.show()

sns.histplot(data=df, x="n_blocks", y="n_categories", bins=20, log_scale=(True, False))
plt.show()

sns.histplot(
    data=df,
    x="n_blocks",
    y="majority_category",
    bins=20,
    hue="singleton",
    log_scale=(True, False),
)
plt.show()

# %%
