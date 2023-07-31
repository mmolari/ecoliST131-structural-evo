# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns

import utils as ut

from Bio import Phylo
from collections import defaultdict

# %%
inf_simple = ut.load_inference(context=False)
inf_context = ut.load_inference(context=True)
tree = ut.load_tree()
bdf = pd.read_csv(ut.expl_fld / "filtered_paths" / "block_stats.csv", index_col=0)
N_strains = len(tree.get_terminals())

# %%


def plot_tree_events(tree, pa_inference, bid, ax):
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
                if et == "gain":
                    node.color = "lime"
                elif et == "loss":
                    node.color = "red"
                else:
                    node.color = "orange"
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
        if pa[nn] == "-":
            return "lightgray"
        else:
            return adj_color[pa[nn]]

    color_tree(tree.root)

    Phylo.draw(
        tree, label_func=label_tree, label_colors=lab_colors, axes=ax, do_show=False
    )
    plt.title(f"block - {bid}")


# %%

svfld = ut.fig_fld / "mugration_trees"
svfld.mkdir(exist_ok=True, parents=True)
for bid, v in inf_simple.items():
    print(bid)

    n = bdf.loc[bid]["count"]
    if (n == 1) or (n == (N_strains - 1)):
        continue

    tree = ut.load_tree()
    pa = v["pa_pattern"]["pa"]
    ev = v["pa_pattern"]["events"]

    fig, axs = plt.subplots(1, 2, figsize=(10, 12))
    plot_tree_events(tree, inf_simple, bid, axs[0])
    plot_tree_events(tree, inf_context, bid, axs[1])
    sns.despine(fig)
    plt.tight_layout()
    plt.savefig(svfld / f"{bid}.png")
    plt.close(fig)
    # break

# %%
