# %%
import numpy as np

import utils as ut

from collections import Counter, defaultdict


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
forbidden_core_flanks = set([e.left.id for e in forbidden_edges]) | set(
    [e.right.id for e in forbidden_edges]
)

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
tree = ut.load_tree()
str_order = [n.name for n in tree.get_terminals()]
# %%

is_dupl = bdf["duplicated"].to_dict()
Ls = bdf["len"].to_dict()

# cmap = plt.get_cmap("rainbow")
# cdict = defaultdict(lambda: cmap(np.random.rand()))
core_anchor = bdf[bdf.core].sort_values("len").index[-1]

fig, ax = plt.subplots(1, 1, figsize=(12, 6))
y = 0
for iso in str_order:
    path = paths[iso]
    x = 0
    B = [node.id for node in path.nodes]
    ca = B.index(core_anchor)
    flip = not path.nodes[ca].strand
    B = np.roll(B, -ca - flip)
    if flip:
        B == B[::-1]
    for bid in B:
        # if Ls[bid] < 500:
        #     continue
        if is_dupl[bid]:
            continue

        l = Ls[bid]
        if bid in forbidden_core_flanks:
            c = "red"
            # c = cdict[bid]
        elif bid in forbidden_blocks:
            c = "pink"
        elif is_core[bid]:
            c = "lightgreen"
        # elif is_dupl[bid]:
        #     c = "white"
        else:
            c = "blue"
        plt.plot([x, x + l], [y, y], color=c)
        x += l
        # if x > 1e5:
        #     break
    y += 1

    # if y > 6:
    #     break
ax.set_xlabel("Genome position (bp)")
ax.set_yticks([])
sns.despine(left=True)
plt.tight_layout()
plt.savefig(ut.fig_fld / "2_synt_breaks_paths.png")
plt.show()

# %%
