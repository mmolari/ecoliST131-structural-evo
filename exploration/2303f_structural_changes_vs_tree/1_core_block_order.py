# %%
import pathlib
import copy

import pypangraph as pp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import Phylo
from collections import defaultdict

fig_fld = pathlib.Path("figs")
fig_fld.mkdir(exist_ok=True)


def svfig(svname):
    for suffix in ["png", "pdf"]:
        plt.savefig(fig_fld / f"{svname}.{suffix}", dpi=300)


# %%


class Block:
    def __init__(self, bid: str, s: bool):
        self.id = bid
        self.s = s

    def __repr__(self):
        sign = "+" if self.s else "-"
        return f"[{self.id}|{sign}]"

    def __eq__(self, b):
        return (self.id == b.id) and (self.s == b.s)

    def __hash__(self):
        return hash((self.id, self.s))

    def invert(self):
        return Block(self.id, not self.s)


class Edge:
    def __init__(self, b1: Block, b2: Block):
        self.b = (b1, b2)

    def __eq__(self, e):
        if self.b == e.b:
            return True
        elif self.invert().b == e.b:
            return True
        else:
            return False

    def __repr__(self):
        return f"({self.b[0]} - {self.b[1]})"

    def __hash__(self):
        e = self.invert()
        return hash((self.b[0], self.b[1])) ^ hash((e.b[0], e.b[1]))

    def invert(self):
        return Edge(self.b[1].invert(), self.b[0].invert())


# %%
pangraph_file = "../../results/ST131/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)

df = pan.to_blockstats_df()
core = df["core"].to_dict()
Ls = df["len"].to_dict()

tree = Phylo.read(
    "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk", format="newick"
)
tree.ladderize()

# %%

sns.histplot(data=df, x="len", hue="core", log_scale=True)
plt.title("block length distribution")
svfig("block_size_distr")
plt.show()

# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 13))
Phylo.draw(
    tree,
    label_func=lambda x: x.name.removeprefix("NZ_")
    if x in tree.get_terminals()
    else "",
    axes=ax,
    do_show=False,
)
for k in ["top", "right", "left"]:
    ax.spines[k].set_visible(False)
ax.set_yticks([])
ax.set_ylabel("")
plt.tight_layout()
svfig("tree")
plt.show()

# %%

paths = {}
for pt in pan.paths:
    bids = pt.block_ids
    bs = pt.block_strands
    paths[pt.name] = [Block(bid, s) for bid, s in zip(bids, bs) if core[bid]]

# longest core block:
Bmax_id = df[df.core].sort_values("len").iloc[-1].name
Bmax = Block(Bmax_id, True)
print(Bmax)

anchor_paths = {}
for k, p in paths.items():
    original_p = paths[k]
    if Bmax in original_p:
        p = copy.copy(original_p)
    else:
        p = [b.invert() for b in original_p[::-1]]
    idx = p.index(Bmax)
    p = np.roll(p, -idx)
    anchor_paths[k] = p
# %%

fig, axs = plt.subplots(
    1, 2, figsize=(12, 8), gridspec_kw={"width_ratios": [1, 4]}, sharey=True
)

ax = axs[0]
Phylo.draw(
    tree,
    label_func=lambda x: "",
    axes=ax,
    do_show=False,
)
# for k in ["top", "right", "left"]:
#     ax.spines[k].set_visible(False)
# ax.set_yticks([])
# ax.set_ylabel("")

ax = axs[1]
cmapf = plt.get_cmap("Blues")
cmapr = plt.get_cmap("Reds")
bcolf = defaultdict(lambda: cmapf(np.random.rand()))
bcolr = defaultdict(lambda: cmapr(np.random.rand()))
bdispl = defaultdict(lambda: np.random.rand() * 2 - 1)

leaves_keys = [l.name for l in tree.get_terminals()]
for n, k in enumerate(leaves_keys):
    p = anchor_paths[k]
    x = 0
    for b in p:
        l = Ls[b.id]
        y = 1 + n + bdispl[b.id] * 0.2
        c = bcolf[b.id] if b.s else bcolr[b.id]
        ax.plot([x, x + l], [y, y], color=c)
        x += l
ax.set_xlabel("core genome alignment")
# ax.set_yticks([])
# for k in ["top", "right", "left"]:
#     ax.spines[k].set_visible(False)
sns.despine()

plt.tight_layout()
svfig("color_paths")
plt.show()
# %%

edge_repr = {}
for k, blks in paths.items():
    edge_repr[k] = [Edge(b1, b2) for b1, b2 in zip(blks, np.roll(blks, -1))]

# %%
edge_strains = defaultdict(set)
for k, edges in edge_repr.items():
    for e in edges:
        edge_strains[e].add(k)

strains = set(pan.strains())

# edges to be merged and edges with a breakpoint
merges = [e for e, es in edge_strains.items() if es == strains]
breakpoints = [e for e, es in edge_strains.items() if es != strains]
# %%

comp_paths = copy.deepcopy(anchor_paths)
Ls = df["len"].to_dict()

merged_bid = {}
for e in merges:
    print(e)
    bid1 = e.b[0].id
    bid2 = e.b[1].id
    s01 = e.b[0].s
    s02 = e.b[1].s
    m1, m2 = bid1 in merged_bid, bid2 in merged_bid
    bid1 = merged_bid[bid1] if m1 else bid1
    bid2 = merged_bid[bid2] if m2 else bid2
    bid = "|".join([bid1, bid2])
    L = Ls[bid1] + Ls[bid2]
    for k in comp_paths:
        p = list(comp_paths[k])

        s1 = Block(bid1, True) in p
        s2 = Block(bid2, True) in p
        s = s1 == s01

        b1 = Block(bid1, s1)
        b2 = Block(bid2, s2)
        idx1 = p.index(b1)
        p.pop(idx1)
        idx2 = p.index(b2)
        p.pop(idx2)
        p.insert(idx2, Block(bid, s))
        comp_paths[k] = p
    Ls[bid] = L
    for x in bid1.split("|") + bid2.split("|"):
        merged_bid[x] = bid


# %%

# %%

fig, axs = plt.subplots(
    1, 2, figsize=(12, 8), gridspec_kw={"width_ratios": [1, 4]}, sharey=True
)

ax = axs[0]
Phylo.draw(
    tree,
    label_func=lambda x: "",
    axes=ax,
    do_show=False,
)
# for k in ["top", "right", "left"]:
#     ax.spines[k].set_visible(False)
# ax.set_yticks([])
# ax.set_ylabel("")

ax = axs[1]
cmap = plt.get_cmap("jet")
bcol = defaultdict(lambda: cmap(np.random.rand()))
bdispl = defaultdict(lambda: np.random.rand() * 2 - 1)

leaves_keys = [l.name for l in tree.get_terminals()]
for n, k in enumerate(leaves_keys):
    p = comp_paths[k]
    x = 0
    for b in p:
        l = Ls[b.id]
        y = 1 + n + bdispl[b.id] * 0.3
        c = bcol[b.id]
        ax.plot([x, x + l], [y, y], color=c)
        if b.s:
            ax.plot([x, x], [y, 1 + n + 0.3], color=c)
        else:
            ax.plot([x + l, x + l], [y, 1 + n - 0.3], color=c)
        x += l
ax.set_xlabel("core genome alignment")
# ax.set_yticks([])
# for k in ["top", "right", "left"]:
#     ax.spines[k].set_visible(False)
sns.despine()

plt.tight_layout()
svfig("color_paths_compressed")
plt.show()
# %%
len(breakpoints)
# %%

for br in breakpoints:
    A = set(edge_strains[br])
    B = set(strains) - A
    if len(A) < len(B):
        print(br, A)
    else:
        print(br, B)
# %%
