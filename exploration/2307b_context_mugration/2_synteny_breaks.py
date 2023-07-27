# %%
import numpy as np

import utils as ut

from collections import Counter


# %% 2a: clean paths

pan = ut.load_pangraph()
bdf = pan.to_blockstats_df()
is_core = bdf["core"].to_dict()
strains = pan.strains()
N = len(strains)

# core-block paths
core_paths = ut.pangraph_to_path_dict(pan)
core_paths = ut.filter_paths(core_paths, lambda bid: is_core[bid])

# %%
# core-edges counter
Edges = []
for iso, path in core_paths.items():
    next_path = np.roll(path.nodes, -1)
    for l, r in zip(path.nodes, next_path):
        edge = ut.Edge(l, r)
        Edges.append(edge)
Edges = Counter(Edges)

# %%
forbidden_edges = set([e for e, n in Edges.items() if n < N])


# %%
paths = ut.pangraph_to_path_dict(pan)
adj_paths = {}
for iso, path in paths.items():
    adj_paths[iso] = ut.to_core_adjacencies(path.nodes, is_core)
# %%

forbidden_blocks = set()
forbidden_edge_ct = Counter()
for iso, Js in adj_paths.items():
    for j in Js:
        bid = j.center.id
        if is_core[bid]:
            print("NEVER HAPPENS")
            continue
        edge = ut.Edge(j.left, j.right)
        if edge in forbidden_edges:
            forbidden_blocks.update([bid])
            forbidden_edge_ct.update([edge])
forbidden_blocks

# %%

bdf["synt. break"] = np.isin(bdf.index, list(forbidden_blocks))
# %%
bdf.groupby(["core", "duplicated"])["synt. break"].value_counts()

# save list of forbidden blocks
np.savetxt(ut.expl_fld / "forbidden_blocks.txt", list(forbidden_blocks), fmt="%s")

# %%

import seaborn as sns
import matplotlib.pyplot as plt

mask = (~bdf["core"]) & (~bdf["duplicated"])
sdf = bdf[mask]

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

ax = axs[0]
sns.histplot(
    data=sdf,
    x="n. strains",
    hue="synt. break",
    ax=ax,
    bins=np.arange(len(strains) + 1) + 0.5,
    element="step",
    # weights="len",
)
ax.set_xlabel("Number of strains")
ax.set_ylabel("Number of accessory non-duplicated blocks")

ax = axs[1]
sns.histplot(
    data=sdf,
    x="len",
    hue="synt. break",
    ax=ax,
    log_scale=True,
    element="step",
)
ax.set_xlabel("Block length (bp)")
ax.set_ylabel("Number of accessory non-duplicated blocks")

plt.tight_layout()
plt.savefig(ut.fig_fld / "2_synt_breaks.png")
plt.show()

# %%
