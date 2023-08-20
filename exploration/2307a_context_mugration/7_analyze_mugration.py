# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import itertools as itt

import utils as ut

from Bio import Phylo
from collections import defaultdict

# %%
inf_context = ut.load_inference(context=True)
tree = ut.load_tree()
bdf = pd.read_csv(ut.expl_fld / "filtered_paths" / "block_stats.csv", index_col=0)
N_strains = len(tree.get_terminals())
leaves = {n.name: n for n in tree.get_terminals()}
leaves_names = list(leaves.keys())

# %%
df = []
for bid, row in bdf[~bdf["core"]].iterrows():
    res = {"bid": bid}

    info = inf_context[bid]
    pa = info["pa_pattern"]["pa"]
    ev = info["pa_pattern"]["events"]

    # number of leaves with/without the block
    vals = np.array([pattern for node, pattern in pa.items() if node in leaves_names])
    res["without"] = (vals == "-").sum()
    res["with"] = (vals != "-").sum()
    res["contexts"] = len(set(vals) - set(["-"]))

    # number of events
    g, l, m = 0, 0, 0
    for node, e in ev:
        if e == "gain":
            g += 1
        elif e == "loss":
            l += 1
        else:
            m += 1
    res["gain"] = g
    res["loss"] = l
    res["move"] = m

    # distances
    ctxts = defaultdict(list)
    for leaf, pattern in pa.items():
        if leaf not in leaves_names:
            continue
        if pattern != "-":
            ctxts[pattern].append(leaf)
    dist_intra = None
    dist_inter = None
    dist_all = None

    # distance all
    all_leaves_with = sum(ctxts.values(), [])
    if len(all_leaves_with) > 1:
        dist_all = 0
        for l1, l2 in itt.combinations(all_leaves_with, 2):
            dist_all += tree.distance(leaves[l1], leaves[l2])
        dist_all /= len(all_leaves_with) * (len(all_leaves_with) - 1) / 2

    if (res["contexts"] > 1) and (res["contexts"] < len(all_leaves_with)):
        dist_intra = 0
        dist_intra_n = 0

        # dist intra
        for ctxt, leaves_with in ctxts.items():
            for l1, l2 in itt.combinations(leaves_with, 2):
                dist_intra += tree.distance(leaves[l1], leaves[l2])
                dist_intra_n += 1
        dist_intra /= dist_intra_n

    if res["contexts"] > 1:
        dist_inter = 0
        dist_inter_n = 0

        # dist inter
        for k1, k2 in itt.combinations(ctxts.keys(), 2):
            for l1 in ctxts[k1]:
                for l2 in ctxts[k2]:
                    dist_inter += tree.distance(leaves[l1], leaves[l2])
                    dist_inter_n += 1
        dist_inter /= dist_inter_n

    res["dist_all"] = dist_all
    res["dist_intra"] = dist_intra
    res["dist_inter"] = dist_inter
    df.append(res)

df = pd.DataFrame(df).set_index("bid")
# %%
df = df.sort_values("with")
df = df.merge(bdf, left_index=True, right_index=True)
df.to_csv(ut.expl_fld / "filtered_paths" / "block_stats_context.csv")

# %%

sns.histplot(data=df, x="with", hue="contexts", multiple="stack")
plt.show()


sns.histplot(data=df, x="len", hue="contexts", multiple="stack", log_scale=True)
plt.show()

# %%
mask = ~df["dist_intra"].isna()
df[mask]
sns.scatterplot(data=df[mask], x="dist_intra", y="dist_inter", hue="with")
plt.plot([0, 5e-5], [0, 5e-5], color="gray")
plt.show()
# %%

mask = df["dist_intra"] > df["dist_inter"]
df[mask].sort_index()

# %%
