# %%

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import utils as ut

from collections import defaultdict


CORE_LEN_THR = 500
SHORT_THR = 200
# %%
pan = ut.load_pangraph()
paths = ut.pangraph_to_path_dict(pan)
bdf = pan.to_blockstats_df()
is_core = bdf["core"].to_dict()
is_dupl = bdf["duplicated"].to_dict()
block_len = bdf["len"].to_dict()
strains = pan.strains()


def core_filter(bid):
    if not is_core[bid]:
        return False
    if block_len[bid] < CORE_LEN_THR:
        return False
    return True


# %%

with open(ut.expl_fld / "backbone_edges.txt", "r") as f:
    backbone_edges = [ut.Edge.from_str_id(l.strip()) for l in f.readlines()]


# %%
def to_core_junctions(path, is_core):
    """Given a path, returns a dictionary of core edges -> core junctions"""
    left_core = None
    right_core = None
    first_core = None
    CJs = {}
    center = []

    # set index to first core node
    i = 0
    while True:
        node = path.nodes[i]
        i += 1
        core = core_filter(node.id)
        if core:
            first_core = node
            left_core = node
            break

    L = len(path.nodes)
    while True:
        node = path.nodes[i]
        core = core_filter(node.id)
        # print(f"i = {i}, node = {node}, core = {core}")
        if core:
            # print(f"left_core = {left_core}, center = {center}, right_core = {node}")
            right_core = node
            J = ut.Junction(left_core, ut.Path(center), right_core)
            E = ut.Edge(left_core, right_core)
            CJs[E] = J
            center = []
            if right_core == first_core:
                break
            else:
                left_core = right_core
        else:
            center.append(node)

        i = (i + 1) % L

    return CJs


CJs = {}
for iso, path in paths.items():
    CJs[iso] = to_core_junctions(path, is_core)
    # break

# %%
# questions:
# - how long on average?
# - what length variance?
# - what is the similarity?
# - how many blocks?

df = {}

for e in backbone_edges:
    df[e] = defaultdict(float)
    for iso in strains:
        cj = CJs[iso][e]

        Bs = [b.id for b in cj.center.nodes]

        df[e]["avg_len"] += sum([block_len[b] for b in Bs])
        df[e]["avg_len_dupl"] += sum([block_len[b] for b in Bs if is_dupl[b]])
        df[e]["avg_n_blocks"] += len(Bs)
        df[e]["avg_n_dupl"] += sum([is_dupl[b] for b in Bs])
        df[e]["avg_n_core"] += sum([is_core[b] for b in Bs])
        df[e]["avg_n_short"] += sum([block_len[b] < 200 for b in Bs])
        df[e]["n_strains"] += 1 * (len(Bs) > 0)
        df[e]["len_max"] = max(df[e]["len_max"], sum([block_len[b] for b in Bs]))
        if not "len_min" in df[e]:
            df[e]["len_min"] = sum([block_len[b] for b in Bs])
        else:
            df[e]["len_min"] = min(df[e]["len_min"], sum([block_len[b] for b in Bs]))

for e in backbone_edges:
    df[e]["avg_len"] /= len(strains)
    df[e]["avg_len_dupl"] /= len(strains)
    df[e]["avg_n_blocks"] /= len(strains)
    df[e]["avg_n_dupl"] /= len(strains)
    df[e]["avg_n_core"] /= len(strains)
    df[e]["avg_n_short"] /= len(strains)

    for iso in strains:
        cj = CJs[iso][e]

        Bs = [b.id for b in cj.center.nodes]

        df[e]["len_var"] += (sum([block_len[b] for b in Bs]) - df[e]["avg_len"]) ** 2
        df[e]["len_dupl_var"] += (
            sum([block_len[b] for b in Bs if is_dupl[b]]) - df[e]["avg_len_dupl"]
        ) ** 2
        df[e]["n_blocks_var"] += (len(Bs) - df[e]["avg_n_blocks"]) ** 2
        df[e]["n_dupl_var"] += (
            sum([is_dupl[b] for b in Bs]) - df[e]["avg_n_dupl"]
        ) ** 2

df = pd.DataFrame(df).T
df = df.sort_values("avg_len", ascending=False)
df
# %%
tree = ut.load_tree()
tree.ladderize()
str_order = [n.name for n in tree.get_terminals()]


def plot_row(ax, nodes, i, colors, lengths):
    x = 0
    for node in nodes:
        color = colors[node.id]
        l = lengths[node.id]
        w = is_dupl[node.id] * 3 + 3
        ax.plot([x, x + l], [i, i], color=color, linewidth=w)

        # dy = 0.3
        # xc = x if node.strand else x + l
        # ax.plot([x, x], [i - dy, i + dy], color=color, linewidth=4)
        # if node.strand:
        #     x_arrow = x + l
        #     arrow_marker = "4"
        # else:
        #     x_arrow = x
        #     arrow_marker = "3"
        # ax.scatter(x_arrow, i, color=color, marker=arrow_marker)

        x += l


def plot_gap(ax, Js, e, str_order, lengths):
    cmap = mpl.cm.get_cmap("rainbow")
    colors = defaultdict(lambda: cmap(np.random.rand()))

    for i, iso in enumerate(str_order):
        j = Js[iso]
        if j.left != e.left:
            j = j.invert()
            assert j.left == e.left
        plot_row(ax, j.center.nodes, i, colors, lengths)

    ax.set_yticks(np.arange(len(str_order)))
    ax.set_yticklabels(str_order)
    ax.set_xlabel("position")


svfld = ut.fig_fld / "core_junctions"
svfld.mkdir(exist_ok=True)

for e in df.index[df["n_strains"] > 1]:
    Js = {iso: CJs[iso][e] for iso in strains}
    fig, ax = plt.subplots(figsize=(20, 20))
    plot_gap(ax, Js, e, str_order, block_len)
    plt.title(f"{e}")
    plt.tight_layout()
    plt.savefig(
        svfld / f"{e.left.to_str_id()}__{e.right.to_str_id()}.png", facecolor="w"
    )
    plt.show()
# %%
save_df = df.copy()
save_df.index = [e.to_str_id() for e in save_df.index]
save_df.to_csv(ut.expl_fld / "core_junctions_df.csv")

# %%
