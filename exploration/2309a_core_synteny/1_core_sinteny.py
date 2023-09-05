# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
from Bio import Phylo
from collections import Counter, defaultdict
import argparse
import pathlib
import utils as ut


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Given a pangraph, generates a list of core edges and their frequencies.
        Only core-blocks longer than the specified threshold are considered.
        """
    )
    parser.add_argument("--pangraph", type=str, required=True)
    parser.add_argument("--tree", type=str, required=True)
    parser.add_argument("--fig_synt", type=str, required=True)
    parser.add_argument("--fig_tree", type=str, required=True)
    parser.add_argument("--fig_blocks", type=str, required=True)
    return parser.parse_args()


fig_fld = pathlib.Path(
    "/home/marco/ownCloud/neherlab/code/pangenome-evo/exploration/2309a_core_synteny/figs"
)
pg_file = "/home/marco/ownCloud/neherlab/code/pangenome-evo/results/ST131/pangraph/asm20-100-5-polished.json"
tree_file = "/home/marco/ownCloud/neherlab/code/pangenome-evo/results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
svfig_synt = fig_fld / "core_synteny.pdf"
svfig_tree = fig_fld / "tree.pdf"
svfig_blocks = fig_fld / "blocks.pdf"

# %%
pan = pp.Pangraph.load_json(pg_file)
bdf = pan.to_blockstats_df()
# %%


def keep_f(bid):
    return bdf.loc[bid, "core"]


paths = ut.pangraph_to_path_dict(pan)
paths = ut.filter_paths(paths, keep_f)
bdf = bdf[bdf["core"]].copy().drop(["core", "duplicated", "n. strains"], axis=1)

# %%


def path_edge_count(paths):
    """Count internal edges of paths"""
    ct = Counter()
    for iso, p in paths.items():
        L = len(p.nodes)
        for i in range(L):
            e = ut.Edge(p.nodes[i], p.nodes[(i + 1) % L])
            ct.update([e])
    return dict(ct)


def find_mergers(edge_ct, block_ct):
    """Create a dictionary source -> sinks of block-ids to be merged"""
    mergers = {}
    for e, ec in edge_ct.items():
        bl, br = e.left.id, e.right.id
        if (ec == block_ct[bl]) and (ec == block_ct[br]):
            # merge
            if bl in mergers:
                mergers[br] = mergers[bl]
            elif br in mergers:
                mergers[bl] = mergers[br]
            else:
                mergers[br] = bl
    return mergers


def perform_mergers(mergers, paths, block_df):
    """Performs selected mergers in paths and block stats"""
    for source, sink in mergers.items():
        # modify block count dataframe
        block_df.loc[sink, "len"] += block_df["len"][source]
        assert block_df["count"][sink] == block_df["count"][source]
        block_df.drop(source, inplace=True)

    removed_blocks = set(mergers.keys()) - set(mergers.values())

    def keep(bid):
        return bid not in removed_blocks

    return ut.filter_paths(paths, keep)


# %%
ec = path_edge_count(paths)
mg = find_mergers(ec, bdf["count"].to_dict())
paths = perform_mergers(mg, paths, bdf)

# %%
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

# %%

path_cats = defaultdict(list)
for iso, path in paths.items():
    path_cats[path].append(iso)

path_cats = dict(path_cats)
print(f"n. different paths = {len(path_cats)}")
abundances = [len(isos) for isos in path_cats.values()]
print(f"abundances = {abundances}")
common_path = max(path_cats, key=lambda p: len(path_cats[p]))
strand_common = {node.id: node.strand for node in common_path.nodes}
# %%

fig, ax = plt.subplots(figsize=(10, 10))

y = 0
# cmap = plt.cm.get_cmap("tab20b")(range(len(bdf)))
cmap = plt.cm.get_cmap("nipy_spectral")(np.linspace(0, 1, len(bdf)))
colors = defaultdict(lambda: cmap[len(colors)])


def cmap_iso_generator():
    cm = plt.cm.get_cmap("tab20")
    i = 0
    while True:
        yield cm(i)
        i += 1


cmap_iso = cmap_iso_generator()


iso_color = {}
xpos = defaultdict(list)
ylabels = []
yticks = []
for path, isolates in path_cats.items():
    for i, node in enumerate(path.nodes):
        xpos[node.id].append(i + 0.5)
        bid = node.id
        c = colors[bid]
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
        yl = ", ".join(isolates)

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
plt.savefig(svfig_synt)
plt.show()

# %%

tree = Phylo.read(tree_file, "newick")

fig, ax = plt.subplots(figsize=(3.5, 10))


Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    show_confidence=False,
    label_func=lambda x: x.name if x.name in iso_color else "",
    label_colors=lambda name: iso_color[name] if name in iso_color else "k",
)
sns.despine()
plt.tight_layout()
plt.savefig(svfig_tree)
plt.show()

# %%

# block lengths

fig, ax = plt.subplots(figsize=(10, 2))

y = 0
yticks, ylabels = [], []
for node in common_path.nodes:
    bid = node.id
    L = bdf.loc[bid, "len"]
    ax.bar(y, L, color=colors[bid], edgecolor="dimgray")
    ylabels.append(f"{int(np.round(L/1000))} kbp")
    yticks.append(y)
    y += 1
ax.set_xticks(yticks)
ax.set_xticklabels(ylabels, rotation=90)
ax.set_ylim(1e3, 1e6)
ax.set_yscale("log")
sns.despine()
ax.grid(which="major", axis="y", alpha=0.6)
# ax.grid(which="minor", axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(svfig_blocks)
plt.show()

# %%

# total tree length
tbl = tree.total_branch_length()
n_events = 15
aln_len = 3772479

print(tbl * aln_len)
# %%
