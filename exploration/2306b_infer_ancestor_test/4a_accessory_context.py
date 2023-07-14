# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pypangraph as pp
import utils4 as ut

from Bio import Phylo
from collections import defaultdict, Counter


# %%

# load pangraph
pan = ut.load_pangraph()
bdf = pan.to_blockstats_df()
is_core = bdf["core"].to_dict()
is_dupl = bdf["duplicated"].to_dict()
bl_Ls = bdf["len"].to_dict()

# load P/A inference
pa_inference = ut.load_pa_inference()

# load tree
tree = ut.load_nodenamed_tree()

fig_fld = ut.fig_fld / "4"
fig_fld.mkdir(exist_ok=True)
# %%

LEN_THRESH = 500


def filter_func(bid):
    if bl_Ls[bid] < LEN_THRESH:
        return False
    if is_dupl[bid]:
        return False
    return True


def to_nodelist(path):
    B = path.block_ids
    S = path.block_strands
    nodes = [ut.Node(b, s) for b, s in zip(B, S)]
    return nodes


def to_acc_adjacencies(nodes):
    adj = {}
    N = len(nodes)

    # appen initial blocks to the end of the list for periodic boundary conditions
    nodes = [nodes[-1]] + list(nodes) + [nodes[0]]

    for i in range(1, N + 1):
        n = nodes[i]
        if is_core[n.id]:
            continue
        adj[n.id] = ut.Junction(nodes[i - 1], n, nodes[i + 1])

    return adj


def to_core_adjacencies(nodes):
    adj = {}
    N = len(nodes)

    for i, n in enumerate(nodes):
        if is_core[n.id]:
            continue
        bef, aft = None, None
        bi, ai = i, i
        while bef is None:
            bi = (bi - 1) % N
            if is_core[nodes[bi].id]:
                bef = nodes[bi]
        while aft is None:
            ai = (ai + 1) % N
            if is_core[nodes[ai].id]:
                aft = nodes[ai]
        adj[n.id] = ut.Junction(bef, n, aft)

    return adj


Adj = {}

for path in pan.paths:
    name = path.name
    print(f"processing {name}")
    nodes = to_nodelist(path)
    print(f"  {len(nodes)} blocks")

    # filter out short and duplicated blocks
    mask = np.array([filter_func(n.id) for n in nodes])
    nodes = np.array(nodes)[mask]
    print(f"  {np.sum(~mask)} blocks filtered")
    print(f"  {len(nodes)} blocks remaining")

    adj_acc = to_acc_adjacencies(nodes)
    adj_core = to_core_adjacencies(nodes)

    Adj[name] = {"acc": adj_acc, "core": adj_core}

    # filter paths

# %%

ctr_a, ctr_c = defaultdict(Counter), defaultdict(Counter)
for iso in pan.strains():
    adj_core, adj_acc = Adj[iso]["core"], Adj[iso]["acc"]
    for k, j in adj_core.items():
        ctr_c[k].update([j])
    for k, j in adj_acc.items():
        ctr_a[k].update([j])

ctr_df = []
for k in ctr_c:
    c_ct = ctr_c[k]
    a_ct = ctr_a[k]
    ctr_df.append(
        {
            "block": k,
            "core_ctx": len(c_ct),
            "acc_ctx": len(a_ct),
            "N": sum(c_ct.values()),
        }
    )
ctr_df = pd.DataFrame(ctr_df)
ctr_df.set_index("block", inplace=True)

# %%

mask = (ctr_df["core_ctx"] > 1) | (ctr_df["acc_ctx"] > 1)
fig, ax = plt.subplots(figsize=(4, 4))
sns.histplot(
    data=ctr_df,
    x="core_ctx",
    y="acc_ctx",
    bins=[np.arange(12) + 0.5] * 2,
    cbar=True,
    ax=ax,
    norm=mpl.colors.LogNorm(),
    vmin=None,
    vmax=None,
    cmap="rainbow",
    cbar_kws={"label": "n. blocks"},
)
ax.set_xlabel("n. core contexts")
ax.set_ylabel("n. accessory contexts")
N_blocks = len(ctr_df)
N_nontrivial = len(ctr_df[mask])
ax.set_title(f"n. blocks: {N_blocks}, n. multiple-contexts: {N_nontrivial}")
plt.tight_layout()
plt.savefig(fig_fld / "core_vs_acc_ctx.png")
plt.show()

# %%
mask = ctr_df["core_ctx"] == ctr_df["acc_ctx"]
ctr_df[mask]

# %%

mask = (ctr_df["core_ctx"] == 2) & (ctr_df["acc_ctx"] == 2)
ctr_df[mask]

# %%


def plot_tree_events(tree, pa_inference, bid, Adj):
    info = pa_inference[bid]
    pa = info["pa_pattern"]["pa"]
    ev = info["pa_pattern"]["events"]

    cm = iter(mpl.cm.tab10.colors)

    def next_color():
        return next(cm)

    adj_color = defaultdict(next_color)

    def color_tree(node):
        for nn, et in ev:
            if node.name == nn:
                node.color = "lime" if et == "gain" else "red"
        if node.color is None:
            node.color = "black"
        for c in node.clades:
            color_tree(c)

    def label_tree(node):
        if node.is_terminal():
            return node.name
        else:
            return ""

    def lab_colors(nn):
        if len(nn) == 0:
            return None
        if pa[nn] == "P":
            j = Adj[nn]["core"][bid]
            return adj_color[j]
        else:
            return "lightgray"

    color_tree(tree.root)

    fig, ax = plt.subplots(1, 1, figsize=(10, 12))
    Phylo.draw(
        tree, label_func=label_tree, label_colors=lab_colors, axes=ax, do_show=False
    )
    plt.title(f"block - {bid}")
    return fig, ax


fig_subfld = fig_fld / "context_trees"
fig_subfld.mkdir(exist_ok=True)

mask = (ctr_df["core_ctx"] > 1) & (ctr_df["acc_ctx"] == ctr_df["core_ctx"])
for bid in ctr_df[mask].index:
    tree = ut.load_nodenamed_tree()
    fig, ax = plot_tree_events(tree, pa_inference, bid, Adj)
    depth = bdf.loc[bid]["count"]
    length = bl_Ls[bid]
    ctxts = ctr_df.loc[bid]["core_ctx"]
    ax.set_title(f"block {bid} | {length} bp | depth = {depth} | contexts = {ctxts}")
    plt.tight_layout()
    plt.savefig(fig_subfld / f"block_{bid}_{depth}_{ctxts}.png")
    plt.show()


# %%
