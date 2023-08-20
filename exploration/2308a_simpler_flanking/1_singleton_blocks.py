# %%

# looking for single accessory blocks always flanked by core blocks

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import utils as ut

from collections import defaultdict
from Bio import Phylo

svfld = ut.fig_fld / "simple_flanking"
svfld.mkdir(exist_ok=True, parents=True)

# %%

pan = ut.load_pangraph()
paths = ut.pangraph_to_path_dict(pan)
bdf = pan.to_blockstats_df()
is_core = bdf["core"]

# %%
Js = {iso: ut.to_junctions(path.nodes) for iso, path in paths.items()}

# %%

# look for non-duplicated accessory blocks
mask = (~bdf["core"]) & (~bdf["duplicated"])
candidates = bdf[mask].index.to_numpy()
print(f"n. accessory candidates = {len(candidates)}")

# %%

# check if they are flanked by core blocks
accepted = list(candidates)

for iso, J in Js.items():
    for j in J:
        bid = j.center.id
        if not bid in accepted:
            continue
        l, r = j.left.id, j.right.id
        if not is_core[l] or not is_core[r]:
            accepted.remove(bid)

# %%

sns.scatterplot(
    data=bdf.loc[accepted],
    x="len",
    y="count",
    alpha=0.5,
)
plt.xscale("log")
plt.xlabel("block length (bp)")
plt.ylabel("block count")
plt.tight_layout()
plt.savefig(svfld / "accepted_block_len_vs_count.png")
plt.show()

# %%

# check how many contexts
contexts = defaultdict(dict)
for iso, J in Js.items():
    for j in J:
        bid = j.center.id
        if not bid in accepted:
            continue
        contexts[bid][iso] = j

for bid, C in contexts.items():
    print(f"bid = {bid}, n. contexts = {len(set(C.values()))}")

# %%

# presence/absence pattern
PA = pan.to_blockcount_df()
PA = PA[accepted]

mask = (PA.sum(axis=0) > 1) & (PA.sum(axis=0) < (len(pan.strains()) - 1))
non_singleton = PA.columns[mask]
# %%

for bid in non_singleton:
    pa = PA[bid]
    tree = ut.load_tree()
    fig, ax = plt.subplots(figsize=(10, 10))
    Phylo.draw(
        tree,
        label_func=lambda x: x.name if x in tree.get_terminals() else None,
        label_colors=lambda x: "C2" if pa[x] else "C3",
        do_show=False,
        show_confidence=False,
        axes=ax,
    )
    plt.title(f"bid = {bid}, L = {bdf.loc[bid, 'len']} bp")
    plt.tight_layout()
    plt.savefig(svfld / f"tree_{bid}.png", facecolor="white")
    plt.show()

# %%

bid = "NPEUXKFGMN"

s = pan.blocks[bid].sequence
# from Bio import SeqIO, Seq
# seqfile = ut.expl_fld / "maps" / f"{bid}.fa"
# with open(seqfile, "w") as f:
#     seq = Seq.Seq(s)
#     rec = SeqIO.SeqRecord(seq, id=bid, description="")
#     SeqIO.write(rec, f, "fasta")
print(s)
# %%

# %%
# save PA pattern
data_fld = ut.expl_fld / "singletons"
data_fld.mkdir(exist_ok=True, parents=True)
PA.index.name = "#name"
PA = PA.applymap(lambda x: "-" if x == 0 else "P")
PA.to_csv(data_fld / "PA.csv")


# %%

# infer mugration
import os


pa_inference = {}

Bs = bdf.loc[non_singleton]
Bs = Bs[Bs["len"] >= 500]
Bs = Bs.index.to_list()
for i, bid in enumerate(Bs):
    print(f"processing {bid} \t - \t {i+1}/{len(Bs)}")
    out_dir = data_fld / bid

    # run treetime mugration
    states = data_fld / "PA.csv"

    cmd = f"""
    treetime mugration \
        --tree {ut.named_tree_file} \
        --states {states} \
        --attribute {bid} \
        --outdir {out_dir} \
    """

    os.system(cmd)

    # output file names
    nexus_file = out_dir / "annotated_tree.nexus"
    rates_file = out_dir / "GTR.txt"

    # parse rates
    rates = ut.parse_gtr(out_dir / "GTR.txt")

    # parse inferred presence/absence pattern
    pa_pattern = ut.parse_nexus(nexus_file)

    os.system(f"rm -r {out_dir}")

    pa_inference[bid] = {
        "pa_pattern": pa_pattern,
        "rates": rates,
    }

# write to json
import json

pa_inference_file = data_fld / "infer_pa.json"
with open(pa_inference_file, "w") as f:
    json.dump(pa_inference, f)

# %%
leaves_names = pan.strains()
df = []
for bid in Bs:
    res = {"bid": bid}

    info = pa_inference[bid]
    pa = info["pa_pattern"]["pa"]
    ev = info["pa_pattern"]["events"]

    # number of leaves with/without the block
    vals = np.array([pattern for node, pattern in pa.items() if node in leaves_names])
    res["without"] = (vals == "-").sum()
    res["with"] = (vals != "-").sum()
    res["contexts"] = len(set(vals) - set(["-"]))

    # number of events
    g, l, m = 0, 0, 0
    for node, e in ev:
        if e == "gain":
            g += 1
        elif e == "loss":
            l += 1
    res["gain"] = g
    res["loss"] = l

    df.append(res)

df = pd.DataFrame(df).set_index("bid")
# %%

import matplotlib as mpl


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


PA = pan.to_blockcount_df()[Bs]
for bid in Bs:
    pa = PA[bid]
    tree = ut.load_tree()
    fig, ax = plt.subplots(figsize=(10, 10))
    plot_tree_events(tree, pa_inference, bid, ax)
    plt.tight_layout()
    plt.savefig(svfld / f"tree_{bid}.png", facecolor="white")
    plt.show()
# %%

x = df["without"]
y = df["gain"] + df["loss"]
plt.scatter(x, y, alpha=0.5)
plt.plot([0, 10], [0, 10], "--", color="gray")
plt.xticks(range(0, 22, 2))
plt.xlabel("number of strains without the block")
plt.ylabel("number of events")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(svfld / "events_vs_without.png", facecolor="white")
plt.show()

# %%
cdf = df.merge(bdf, left_index=True, right_index=True)
print(
    cdf[["without", "count", "len", "gain", "loss"]]
    .sort_values("len", ascending=False)
    .to_markdown()
)

# %%
