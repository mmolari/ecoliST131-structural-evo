# %%

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import pypangraph as pp
import utils as ut

from Bio import Phylo
from collections import Counter, defaultdict

e = "KIHCQKPRTH_f__SIKBWUTDMT_f"
pan_file = f"../../results/ST131/backbone_joints/asm20-100-5/joints_pangraph/{e}.json"
tree_file = "data/named_tree.nwk"
filter_len = None
edge = pan_file.split("/")[-1].split(".")[0]

# %%


pan = pp.Pangraph.load_json(pan_file)
bdf = pan.to_blockstats_df()
block_len = bdf["len"].to_dict()
is_core = bdf["core"].to_dict()
is_dupl = bdf["duplicated"].to_dict()

paths = ut.pangraph_to_path_dict(pan)

# optional: clean up paths

if filter_len is not None:
    paths = ut.filter_paths(paths, lambda x: block_len[x] >= filter_len)


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


path_cat = path_categories(paths)

# %%


def plot_row(ax, nodes, y, colors, block_df):
    x = 0
    for node in nodes:
        color = colors[node.id]
        l = block_df["len"][node.id]
        h = block_df["duplicated"][node.id] * 0.2 + 0.4
        core = block_df["core"][node.id]
        edge_color = "black" if node.strand else "red"
        ax.barh(
            y,
            l,
            height=h,
            left=x,
            color=color,
            edgecolor=edge_color,
            hatch="." if core else None,
        )
        x += l


def plot_categories(path_categories, block_df, tree_file):
    # load tree
    tree = Phylo.read(tree_file, "newick")
    tree.ladderize()
    leaves = [n for n in tree.get_terminals()]

    # assign colors to leaves
    C = len(path_categories)

    if C <= 10:
        path_colors = mpl.cm.get_cmap("tab10")(np.arange(C))
    elif C <= 20:
        path_colors = mpl.cm.get_cmap("tab20")(np.arange(C))
    else:
        path_colors = mpl.cm.get_cmap("jet")(np.linspace(0, 1, C))

    strain_color = defaultdict(lambda: "white")
    for i, (_, _, isolates) in enumerate(path_categories):
        for iso in isolates:
            strain_color[iso] = path_colors[i]

    # assign color to blocks
    N_blocks = len(block_df)
    colors = mpl.cm.get_cmap("rainbow")(np.linspace(0, 1, N_blocks))
    np.random.shuffle(colors)
    block_color = {block: colors[i] for i, block in enumerate(block_df.index)}

    fig, axs = plt.subplots(
        1, 2, figsize=(20, 15), gridspec_kw={"width_ratios": [1, 3]}
    )

    ax = axs[0]
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda x: x.name if x in leaves else "",
        label_colors=lambda x: strain_color[x],
    )
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax = axs[1]
    for i, (count, nodes, isolates) in enumerate(path_categories):
        plot_row(ax, nodes, -i, block_color, block_df)
        ax.text(
            0,
            -i + 0.45,
            f"path {i+1} | n = {count}",
            color=path_colors[i],
            fontsize=20,
        )

    ax.set_yticks([])
    ax.set_ylim(-len(path_categories) + 0.5, 0.5)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("position")

    plt.tight_layout()
    return fig, ax


fig, ax = plot_categories(path_cat, bdf, tree_file)
ax.set_title(f"{edge}")

# plt.savefig(
#     svfld / f"{e.left.to_str_id()}__{e.right.to_str_id()}.png", facecolor="w"
# )
# plt.close(fig)
plt.show()

# %%
