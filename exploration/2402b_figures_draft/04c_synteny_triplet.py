# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pypangraph as pp
from Bio import SeqIO, Phylo
import pandas as pd
import pathlib
from collections import defaultdict

fig_fld = pathlib.Path("figs/f04c")
fig_fld.mkdir(parents=True, exist_ok=True)

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(fname, "newick")
tree.ladderize()

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(fname)
bdf = pan.to_blockstats_df()
# %%
fname = "../../results/ST131_ABC/pangraph/genome_lengths.csv"
Ls = pd.read_csv(fname, index_col=0)["length"].to_dict()

fname = "../../exploration/2402b_figures_draft/figs/f04/ST131_ABC/block_stats.csv"
bc = pd.read_csv(fname, index_col=0)
fname = "../../exploration/2402b_figures_draft/figs/f04/ST131_ABC/mergers.csv"
mg = pd.read_csv(fname, index_col=0)
mg = mg["0"].to_dict()

# %%
clade = tree.root[-1][-1][-1][-1][-1][-1][-1][-1][-1][-1][-1][8]
fig, ax = plt.subplots(1, 1, figsize=(2.5, 4))
Phylo.draw(
    clade,
    axes=ax,
    do_show=False,
    label_func=lambda x: x.name if x in clade.get_terminals() else "",
)
ax.set_title("")
ax.set_yticks([])
ax.set_ylabel("")
ax.set_xlabel("")
ax.spines["left"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(fig_fld / "sub_tree.svg")
plt.show()

leaves = [l.name for l in clade.get_terminals()]
# %%
# QQSILILDBT 36253	#7b008c
# ELYKQDZNBM 827	#860097
# KGJXRPELGX 1919	#3800a3


center_block = "ELYKQDZNBM"
thr_dist = 1e5


def circular_distance(x, y, L):
    return np.minimum((x - y) % L, (y - x) % L)


def center_position(block, Bs, Ps):
    idx = list(Bs).index(block)
    p = Ps[idx]
    return p


def pos_dict(Bs, Ss, Ps):
    Pd = defaultdict(list)
    for i, (b, s) in enumerate(zip(Bs, Ss)):
        pos = (Ps[i], Ps[i + 1]) if s else (Ps[i + 1], Ps[i])
        Pd[b].append(pos)
    return Pd


def extract_blocks(path):
    Bs, Ss, Ps = path.block_ids, path.block_strands, path.block_positions
    return Bs, Ss, Ps


gray_cm = mpl.colormaps.get_cmap("gray")
new_colors = defaultdict(lambda: gray_cm(np.random.rand()))
center_pos = {}
fig, axs = plt.subplots(3, 3, figsize=(12, 12), sharex="col", sharey="row")
# equal size for all subplots
# for i in range(3):
#     for j in range(3):
#         axs[i, j].set_aspect("equal")
for i, iso_i in enumerate(leaves):
    for j, iso_j in enumerate(leaves):
        if i <= j:
            if (i - 1 >= 0) and (j <= 2):
                ax = axs[i - 1, j]
                ax.axis("off")
            continue
        ax = axs[i - 1, j]
        path_i = pan.paths[iso_i]
        path_j = pan.paths[iso_j]
        Bs_i, Ss_i, Ps_i = extract_blocks(path_i)
        Bs_j, Ss_j, Ps_j = extract_blocks(path_j)
        Pd_i = pos_dict(Bs_i, Ss_i, Ps_i)
        Pd_j = pos_dict(Bs_j, Ss_j, Ps_j)
        ci = center_position(center_block, Bs_i, Ps_i)
        cj = center_position(center_block, Bs_j, Ps_j)
        center_pos[iso_i] = ci
        center_pos[iso_j] = cj
        for bi, ti in Pd_i.items():
            if bi not in Pd_j:
                continue
            for pi in ti:
                if circular_distance(pi[0], ci, Ls[iso_i]) > thr_dist:
                    continue
                if bi in mg:
                    c = bc["color"][mg[bi]]
                elif bi in bc.index:
                    c = bc["color"][bi]
                else:
                    # c = new_colors[bi]
                    c = "lightgray"

                for pj in Pd_j[bi]:
                    if circular_distance(pj[0], cj, Ls[iso_j]) > thr_dist:
                        continue
                    ax.plot(
                        pj,
                        pi,
                        color=c,
                    )
        if i == 3:
            ax.set_xlabel(iso_j)
        if j == 0:
            ax.set_ylabel(iso_i)
plt.tight_layout()
plt.savefig(fig_fld / "triplet.svg")
plt.show()

# %%


def load_annotations(iso):
    fnames = [
        f"../../results/ST131_ABC/annotations/loc/{tp}.csv"
        for tp in ["defensefinder", "genomad", "integronfinder", "ISEScan"]
    ]
    dfs = []
    for fname in fnames:
        df = pd.read_csv(fname, index_col=0)
        mask = df["iso"] == iso
        df = df[mask]
        dfs.append(df)

    return pd.concat(dfs)


iso = "NZ_CP115324.1"
center = 1e6
width = 1e4
df = load_annotations(iso)

fig, ax = plt.subplots(1, 1, figsize=(8, 3))
y = 1
with open(f"../../data/gbk/{iso}.gbk", "r") as f:
    for record in SeqIO.parse(f, "genbank"):
        for feat in record.features:
            if feat.type == "CDS":
                start = feat.location.start
                if circular_distance(start, center, Ls[iso]) > width:
                    continue
                end = feat.location.end
                color = "gray"
                if "gene" in feat.qualifiers:
                    gene = feat.qualifiers["gene"][0]
                if "product" in feat.qualifiers:
                    product = feat.qualifiers["product"][0]

                y = (y - 0.1) % 1
                ax.plot([start, end], [y, y], "|-", color=color)
                ax.text(
                    start,
                    y,
                    # gene,
                    product,
                    rotation=0,
                    verticalalignment="bottom",
                    fontsize=6,
                )

for _, row in df.iterrows():
    beg, end = row["beg"], row["end"]
    if circular_distance(beg, center, Ls[iso]) > width:
        continue
    tp = row["type"]
    match tp:
        case "IS":
            color = "blue"
        case "CALIN":
            color = "green"
        case _:
            color = "red"
    y = 1.2
    ax.plot([beg, end], [y, y], "|-", color=color)
    ax.text(
        beg,
        y,
        tp,
        rotation=0,
        verticalalignment="bottom",
        fontsize=6,
    )
ax.set_yticks([])
ax.set_ylabel("")
ax.set_xlabel(f"{iso} genome (bp)")
ax.spines["left"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(fig_fld / f"annotations_{iso}.svg")
plt.show()
# %%
