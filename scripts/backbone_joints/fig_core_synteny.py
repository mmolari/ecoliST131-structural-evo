import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pypangraph as pp
from Bio import Phylo
from collections import defaultdict
import argparse
import pathlib
import utils as ut


def armytage_cmap():
    colors = [
        "#f0a3ff",
        "#0075dc",
        "#993f00",
        # "#191919",
        "#4c005c",
        "#005c31",
        "#2bce48",
        "#ffcc99",
        "#808080",
        "#94ffb5",
        "#8f7c00",
        "#9dcc00",
        "#c20088",
        "#003380",
        "#ffa405",
        "#ffa8bb",
        "#426600",
        "#ff0010",
        "#5ef1f2",
        "#00998f",
        # "#e0ff66",
        "#740aff",
        "#990000",
        "#ff5005",
        "#ffff00",
    ]
    for c in colors:
        yield c
    while True:
        yield "k"


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Given a pangraph, generates a list of core edges and their frequencies.
        Only core-blocks longer than the specified threshold are considered.
        """
    )
    parser.add_argument("--pangraph", type=str, required=True)
    parser.add_argument("--tree", type=str, required=True)
    parser.add_argument("--len_thr", type=int, required=True)
    parser.add_argument("--fig_fld", type=str, required=True)
    return parser.parse_args()


def perform_mergers(mergers, paths, block_df):
    """Performs selected mergers in paths and block stats"""
    sink_blocks = set(mergers.values())

    for source, sink in mergers.items():
        # modify block count dataframe
        if source in sink_blocks:
            continue
        # assert source not in kept_blocks, f"{source} in kept_blocks for sink {sink}"
        block_df.loc[sink, "len"] += block_df["len"][source]
        assert block_df["count"][sink] == block_df["count"][source]
        block_df.drop(source, inplace=True)

    discarded_blocks = set(mergers.keys()) - sink_blocks

    def keep(bid):
        return bid not in discarded_blocks

    return ut.filter_paths(paths, keep)


def roll_to_longest_block(paths, bdf):
    longest_bid = bdf["len"].sort_values().index[-1]

    def roll_to_anchor(paths, anchor):
        """Rolls the path so that the anchor is the first node"""
        anchor_node = ut.Node(anchor, True)
        for iso, path in paths.items():
            if not (anchor_node in path.nodes):
                path = path.invert()
            assert anchor_node in path.nodes

            idx = path.nodes.index(anchor_node)
            path.nodes = path.nodes[idx:] + path.nodes[:idx]
            paths[iso] = path

    roll_to_anchor(paths, longest_bid)


def fig_syntey(path_cats, common_path, strand_common, bdf, svname):
    fig, ax = plt.subplots(figsize=(10, 10))

    # cmap = plt.cm.get_cmap("tab20b")(range(len(bdf)))
    cmap = mpl.cm.get_cmap("nipy_spectral")(np.linspace(0, 1, len(bdf)))
    block_colors = defaultdict(lambda: cmap[len(block_colors)])

    cmap_iso = armytage_cmap()

    iso_color = {}
    xpos = defaultdict(list)
    y = 0
    for path, isolates in path_cats.items():
        for i, node in enumerate(path.nodes):
            xpos[node.id].append(i + 0.5)
            bid = node.id
            c = block_colors[bid]
            st = node.strand if strand_common[bid] else not node.strand
            ax.barh(
                y,
                1,
                left=i,
                color=c,
                label=bid,
                edgecolor="dimgray",
                hatch=None if st else "///",
                height=0.3 + 0.5 * (bdf.loc[bid, "len"] > 100000),
            )

        if len(isolates) > 2:
            yl = f"n = {len(isolates)}"
        else:
            yl = "\n".join(isolates)

        if len(isolates) > 10:
            ic = "k"
        else:
            ic = next(cmap_iso)
            for iso in isolates:
                iso_color[iso] = ic
        ax.annotate(
            yl,
            xy=(0, y),
            xytext=(-10, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            color=ic,
        )

        y -= 1

    for bid, X in xpos.items():
        ax.plot(X, -np.arange(len(X)), "k", lw=0.5, zorder=-1)

    # ax.set_yticks(yticks)
    # ax.set_yticklabels(ylabels)
    ax.set_yticks([])
    ax.set_xlim(-10, len(common_path.nodes) + 1)

    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.savefig(str(svname) + ".pdf")
    plt.savefig(str(svname) + ".svg")
    plt.close(fig)
    return iso_color, block_colors


def fig_tree(tree, iso_color, svname):
    fig, ax = plt.subplots(figsize=(3.8, 10))

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda x: "",
        # label_func=lambda x: x.name if x.name in iso_color else "",
        # label_colors=lambda name: iso_color[name] if name in iso_color else "k",
    )

    y_strain = {l.name: n + 1 for n, l in enumerate(tree.get_terminals())}

    xl = ax.get_xlim()[1]
    for iso in iso_color:
        x = tree.distance(iso)
        y = y_strain[iso]
        xf = x + 0.2 * xl
        ax.scatter(xf, y, s=30, c=iso_color[iso], edgecolor="k", zorder=10)
        ax.plot([x, xf], [y, y], c="gray", ls=":")

    sns.despine()
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([0, 5e-5, 1e-4])
    plt.tight_layout()
    plt.savefig(str(svname) + ".pdf")
    plt.savefig(str(svname) + ".svg")
    plt.close()


def fig_blocks(common_path, bdf, block_colors, svname):
    fig, ax = plt.subplots(figsize=(10, 2))

    y = 0
    yticks, ylabels = [], []
    for node in common_path.nodes:
        bid = node.id
        L = bdf.loc[bid, "len"]
        ax.bar(y, L, color=block_colors[bid], edgecolor="dimgray")
        ylabels.append(f"{int(np.round(L/1000))} kbp")
        yticks.append(y)
        y += 1
    ax.set_xticks(yticks)
    ax.set_xticklabels(ylabels, rotation=90)
    ax.set_ylim(500, 1e6)
    ax.set_yscale("log")
    sns.despine()
    ax.grid(which="major", axis="y", alpha=0.6)
    # ax.grid(which="minor", axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(str(svname) + ".pdf")
    plt.savefig(str(svname) + ".svg")
    plt.close()


if __name__ == "__main__":
    args = parse_args()

    # make folder and figure names
    fig_fld = pathlib.Path(args.fig_fld)
    fig_fld.mkdir(exist_ok=True, parents=True)

    svfig_synt = fig_fld / f"core_synteny"
    svfig_tree = fig_fld / f"tree"
    svfig_blocks = fig_fld / f"blocks"

    len_thr = args.len_thr

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)
    bdf = pan.to_blockstats_df()

    # load tree
    tree = Phylo.read(args.tree, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    strains = [l.name for l in tree.get_terminals()]

    # extract paths and filter to core paths
    def keep_f(bid):
        return bdf.loc[bid, "core"] and (bdf.loc[bid, "len"] > len_thr)

    paths = ut.pangraph_to_path_dict(pan)
    paths = ut.filter_paths(paths, keep_f)
    mask = bdf["core"] & (bdf["len"] > len_thr)
    bdf = bdf[mask].copy().drop(["core", "duplicated", "n. strains"], axis=1)

    # merge transitive edges
    mg = ut.find_mergers(paths)
    paths = perform_mergers(mg, paths, bdf)

    # roll to longest block
    roll_to_longest_block(paths, bdf)

    # path categories
    path_cats = defaultdict(list)
    for iso in strains:
        path = paths[iso]
        path_cats[path].append(iso)

    path_cats = dict(path_cats)
    print(f"n. different paths = {len(path_cats)}")
    abundances = [len(isos) for isos in path_cats.values()]
    print(f"abundances = {abundances}")
    # find most common path
    common_path = max(path_cats, key=lambda p: len(path_cats[p]))
    strand_common = {node.id: node.strand for node in common_path.nodes}

    # synteny figure
    iso_color, block_colors = fig_syntey(
        path_cats, common_path, strand_common, bdf, svfig_synt
    )

    # tree figure
    fig_tree(tree, iso_color, svfig_tree)

    # block lengths figure
    fig_blocks(common_path, bdf, block_colors, svfig_blocks)
