# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
from Bio import Phylo
from collections import Counter, defaultdict
import pathlib
import utils as ut


fig_fld = pathlib.Path("figs")
pg_file = f"../../results/ST131/pangraph/asm20-100-5-polished.json"
tree_file = f"../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
svfig_name = fig_fld / f"core_synteny_len.svg"

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
                if br in mergers:
                    source = mergers[br]
                    sink = mergers[bl]
                    for k in mergers:
                        if mergers[k] == source:
                            mergers[k] = sink
                else:
                    mergers[br] = mergers[bl]
            elif br in mergers:
                mergers[bl] = mergers[br]
            else:
                mergers[br] = bl
                mergers[bl] = bl
    return mergers


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

fig, axs = plt.subplots(
    1, 2, figsize=(12, 10), gridspec_kw={"width_ratios": [0.2, 0.8]}
)

ax = axs[1]
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
for path, isolates in sorted(path_cats.items(), key=lambda x: -len(x[1])):
    x = 0
    for i, node in enumerate(path.nodes):
        bid = node.id
        c = colors[bid]
        st = node.strand if strand_common[bid] else not node.strand
        l = bdf.loc[bid, "len"]
        ax.barh(
            y + 0.15 * (st - 0.5),
            l,
            left=x,
            color=c,
            label=bid,
            edgecolor="dimgray",
            hatch=None if st else "///",
            # height=0.3 + 0.5 * (bdf.loc[bid, "len"] > 100000),
            height=0.5,
        )
        xpos[node.id].append(x + l * 0.5)
        xpos[node.id].append(x + l * 0.5)
        x += l

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
    ax.plot(X, -np.arange(len(X)) / 2 + 0.25, color=colors[bid], lw=2, zorder=-1)

# ax.set_yticks(yticks)
# ax.set_yticklabels(ylabels)
ax.set_yticks([])
ax.set_xlabel("core genome (bp)")
# ax.set_xlim(-10, len(common_path.nodes) + 1)

# ax.set_xticks([])
for k in ["top", "right", "left"]:
    ax.spines[k].set_visible(False)


tree = Phylo.read(tree_file, "newick")

ax = axs[0]
Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    show_confidence=False,
    label_func=lambda x: x.name if x.name in iso_color else "",
    label_colors=lambda name: iso_color[name] if name in iso_color else "k",
)
for k in ["top", "right"]:
    ax.spines[k].set_visible(False)
plt.tight_layout()
plt.savefig(svfig_name)
plt.show()

# %%
