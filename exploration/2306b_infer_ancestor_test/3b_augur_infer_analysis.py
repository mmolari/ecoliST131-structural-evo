# %%

import pathlib
import os
import json

import itertools as itt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import pypangraph as pp
import utils3 as ut

from Bio import Phylo
from collections import Counter

fld = pathlib.Path("data/mugration")
fld.mkdir(exist_ok=True)

fig_fld = pathlib.Path("figs/mugration")
fig_fld.mkdir(exist_ok=True)

# %%

prefix = "../../results/ST131/pangraph"
original_tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
tree_file = fld / "tree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

# load tree and pangraph
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()

# %%

# select non-duplicated accessory genes
mask = ~bdf["core"]
mask &= ~bdf["duplicated"]

# list of accessory blocks ids
B_acc = bdf[mask].index.to_list()

# load gene presence/absence inference
pa_inference_file = fld / "infer_pa.json"
with open(pa_inference_file, "r") as f:
    pa_inference = json.load(f)


# %%

# create summary table
summary = []
for b in B_acc:
    info = pa_inference[b]
    pa = info["pa_pattern"]["pa"]
    r = info["rates"]
    ev = info["pa_pattern"]["events"]
    # count event type
    n_gain = len([e for e in ev if e[1] == "gain"])
    n_loss = len([e for e in ev if e[1] == "loss"])
    # count n. of leaves
    pa_count = [v for k, v in pa.items() if not k.startswith("node_")]
    pa_count = Counter(pa_count)
    summary.append(
        {
            "bid": b,
            "n_events": len(ev),
            "n_gains": n_gain,
            "n_losses": n_loss,
            "n_leaves_P": pa_count["P"],
            "n_leaves_A": pa_count["A"],
        }
    )
df = pd.DataFrame(summary).set_index("bid")
df["max_events"] = df[["n_leaves_P", "n_leaves_A"]].min(axis=1)
Ls = bdf["len"].to_dict()
df["len"] = df.apply(lambda x: Ls[x.name], axis=1)


# define category
df["category"] = "multiple-events"
mask = df["n_events"] == 1
mask &= df["max_events"] > 1
df.loc[mask, "category"] = "perfect-support"
mask = df["max_events"] == 1
df.loc[mask, "category"] = "trivial"


df.to_csv(fld / "infer_pa_summary.csv")

# %%
df[["n_events", "n_gains", "n_losses", "max_events"]].value_counts()

# %%

fig, axs = plt.subplots(2, 2, figsize=(10, 8))

ax = axs[0, 0]
bins = np.arange(df["n_leaves_P"].max() + 2) - 0.5
sns.histplot(df, x="n_leaves_P", bins=bins, ax=ax, element="step")
ax.set_xlabel("n. isolates")
ax.set_ylabel("n. blocks")

ax = axs[0, 1]
bins = np.logspace(np.log10(df["len"].min()) - 0.1, np.log10(df["len"].max()) + 0.1, 40)
sns.histplot(df, x="len", bins=bins, ax=ax, element="step")
ax.set_xscale("log")
ax.set_xlabel("block length")
ax.set_ylabel("n. blocks")


ax = axs[1, 0]
M = max(df["n_events"].max(), df["max_events"].max())
bins = np.arange(M + 2) - 0.5
bins = (bins, bins)
sns.histplot(
    df,
    x="n_events",
    y="max_events",
    bins=bins,
    cbar=True,
    cbar_kws={
        # "shrink": 0.5,
        "label": "n. blocks",
    },
    norm=mpl.colors.LogNorm(),
    vmin=None,
    vmax=None,
    # weights="len",
    ax=ax,
)
ax.set_xlabel("n. events")
ax.set_ylabel("n. leaves with minority pattern")
ax.plot([0, M], [0, M], "k--")

ax = axs[1, 1]
M = max(df["n_gains"].max(), df["n_losses"].max())
bins = np.arange(M + 2) - 0.5
bins = (bins, bins)
sns.histplot(
    df,
    x="n_gains",
    y="n_losses",
    bins=bins,
    cbar=True,
    cbar_kws={
        # "shrink": 0.5,
        "label": "n. blocks",
    },
    norm=mpl.colors.LogNorm(),
    vmin=None,
    vmax=None,
    # weights="len",
    ax=ax,
)
ax.set_xlabel("n. gains")
ax.set_ylabel("n. losses")
ax.plot([0, M], [0, M], "k--")
plt.tight_layout()
plt.savefig(fig_fld / "n_events.png")
plt.show()
# %%

mask = df["max_events"] < df["n_events"]
df[mask]

# %%
mask = df["max_events"] < df["n_events"]
for bid in df[mask].index.to_list():
    fig, ax = ut.plot_tree_events(tree_file, pa_inference, bid)
    plt.tight_layout()
    plt.savefig(fig_fld / f"tree_{bid}.png")
    plt.show()


# %%
# tyrosine-type recombinase/integrase, ~ 1kb, almost exact overlap
pan.blocks["WNAPVOJCVV"].sequence

# DNA transfer protein
pan.blocks["JGMLMOYQKZ"].sequence

# tRNA-Ser
pan.blocks["MANREFTSIL"].sequence

# %%

mask = df["max_events"] < df["n_events"]
for bid in df[mask].index.to_list():
    s = pan.blocks[bid].sequence
    qry_fname = f"blast/{bid}.fa"
    with open(qry_fname, "w") as f:
        f.write(">ARODSNUKVQ\n")
        f.write(s)

    cmd = f"""blastn -query {qry_fname} -db blastdb/genomes_db -out blast/{bid}.txt \
        # -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\"
    """
    os.system(cmd)


# %%

# TODO:
# - document this
# - look for branches with good support (1 event / non-trivial)
# - what is the rest?
# - hist of accessory blocks PA frequencies

# %%

# %%

fig, axs = plt.subplots(1, 2, figsize=(10, 4))

ax = axs[0]
bins = np.arange(df["n_leaves_P"].max() + 2) - 0.5
sns.histplot(
    df,
    x="n_leaves_P",
    bins=bins,
    ax=ax,
    hue="category",
    element="step",
    multiple="stack",
)
ax.set_xlabel("n. isolates")
ax.set_ylabel("n. blocks")

ax = axs[1]
bins = np.logspace(np.log10(df["len"].min()) - 0.1, np.log10(df["len"].max()) + 0.1, 40)
sns.histplot(
    df,
    x="len",
    bins=bins,
    ax=ax,
    hue="category",
    element="step",
    multiple="stack",
)
ax.set_xscale("log")
ax.set_xlabel("block length")
ax.set_ylabel("n. blocks")

plt.tight_layout()
plt.savefig(fig_fld / "category_hist.png")
plt.show()

# %%

df.groupby("category")["len"].sum() / df["len"].sum()

# %%
df.groupby("category")["category"].count()

# %%

df.groupby("category")["category"].count() / df.shape[0]

# %%

df[df["category"] == "perfect-support"].sort_values("len", ascending=False)
# %%
bid = "QAICACBTFF"
ut.plot_tree_events(tree_file, pa_inference, bid)
plt.show()
# %%
mask = df["category"] == "perfect-support"
sdf = df[mask]

support = []

for bid in sdf.index.to_list():
    info = pa_inference[bid]
    ev = info["pa_pattern"]["events"]
    support.append(ev[0][0])

support = Counter(support)

tree = Phylo.read(tree_file, "newick")

cmap = plt.get_cmap("cool_r")
norm = mpl.colors.Normalize(vmin=1, vmax=support.most_common()[0][1])


def color_node(node):
    if node.name in support:
        color = cmap(norm(support[node.name]))
        node.color = mpl.colors.to_hex(color)
    else:
        node.color = "gray"
    if not node.is_terminal():
        for c in node.clades:
            color_node(c)


color_node(tree.root)

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: None)
plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="n. events")
plt.tight_layout()
plt.savefig(fig_fld / "tree_support.png")
plt.show()


# %%
