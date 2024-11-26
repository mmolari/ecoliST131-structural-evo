# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
from utils import name_tree_nodes
from Bio import Phylo, SeqIO
import json


res_fld = pathlib.Path("res")

# load gene counts
GC = pd.read_csv(res_fld / "gene_counts.csv", index_col=0)


# load tree:
pangraph_tree_file = (
    "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
)
tree = Phylo.read(pangraph_tree_file, "newick")
name_tree_nodes(tree)
strains = [l.name for l in tree.get_terminals()]

# load gene clusters info
gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"

with open(gene_cluster_json, "r") as f:
    gdf = pd.DataFrame(json.load(f))
gdf["divers"] = gdf["divers"].astype(float)
gdf["count"] = gdf["count"].astype(int)
gdf["geneLen"] = gdf["geneLen"].astype(int)
gdf["event"] = gdf["event"].astype(int)
gdf["core"] = (gdf["count"] == len(strains)) & (gdf["dupli"] == "no")
gdf.set_index("geneId", inplace=True)

mask_core = (GC == GC.iloc[0]).all(axis=0)
GC_acc = GC.loc[:, ~mask_core]

E = pd.read_csv(res_fld / "mugr_events.csv")
E

# %%


def extra_n_dupl(txt):
    extra = 0
    T = txt.strip()
    if len(T) == 0:
        return 0
    for x in T.split("@"):
        print(f"x: {x}")
        k = int(x.split("#")[1])
        extra += k - 1
    return extra


def extract_event_n(df, gdf=gdf):
    is_gain = df["type"] == "gain"
    n_gains = df[is_gain].groupby("gcl").apply(lambda x: x["delta"].sum())
    n_losses = df[~is_gain].groupby("gcl").apply(lambda x: -x["delta"].sum())

    res_df = pd.DataFrame({"n_gains": n_gains, "n_losses": n_losses})
    res_df.fillna(0, inplace=True)
    res_df = res_df.astype(int)

    res_df["n_events"] = res_df["n_gains"] + res_df["n_losses"]

    res_df["n_strains"] = gdf.loc[res_df.index, "count"]
    res_df["n_occ"] = gdf.loc[res_df.index, "dup_detail"].apply(
        lambda x: extra_n_dupl(x)
    )
    res_df["n_occ"] += gdf.loc[res_df.index, "count"]
    res_df["dupli"] = gdf.loc[res_df.index, "dupli"]
    return res_df


# mask = E["branch"].isin(["int_node_1", "int_node_35"])
edf = extract_event_n(E)
edf
# %%

fig_fld = pathlib.Path("figs/f03")
fig_fld.mkdir(exist_ok=True, parents=True)


fig, ax = plt.subplots(figsize=(10, 10))
sns.histplot(
    data=edf,
    y="n_events",
    discrete=True,
    x="n_strains",
    hue="dupli",
    hue_order=["no", "yes"],
    ax=ax,
)
ax.plot([0, 111, 222], [0, 111, 0], "--", color="gray", alpha=0.3)
plt.tight_layout()
# plt.savefig(fig_fld / "n_events_vs_n_strains.png")
plt.show()

# %%
top = edf.sort_values("n_events", ascending=False).head(5)

print(gdf.loc[top.index]["ann"].to_markdown())
print(top.to_markdown())

edf["gains - losses"] = edf["n_gains"] - edf["n_losses"]
dupl_mask = edf["dupli"] == "yes"

fig, axs = plt.subplots(2, 1, figsize=(8, 6))

ax = axs[0]
sns.histplot(edf[~dupl_mask], x="gains - losses", discrete=True, ax=ax, element="step")
ax.axvline(0, color="k", linestyle="--")
ax.set_title("non-duplicated genes")
ax.set_ylabel("n. gene clusters")
ax.set_xlabel("n. gains - n. losses")
ax.set_xlim(-20, 50)

ax = axs[1]
sns.histplot(edf[dupl_mask], x="gains - losses", discrete=True, ax=ax, element="step")
ax.axvline(0, color="k", linestyle="--")
ax.set_title("Duplicated genes")
ax.set_ylabel("n. gene clusters")
ax.set_xlabel("n. gains - n. losses")
ax.set_xlim(-20, 50)

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "gains_minus_losses.png", dpi=300)
plt.savefig(fig_fld / "gains_minus_losses.pdf")
plt.show()


# %%
fig, axs = plt.subplots(1, 2, figsize=(10, 4))

dupl_mask = edf["dupli"] == "yes"

ax = axs[0]
sns.histplot(
    data=edf[~dupl_mask],
    y="n_gains",
    discrete=True,
    x="n_strains",
    ax=ax,
)
ax.plot([0, 111, 222], [0, 111, 0], "--", color="gray", alpha=0.3)
ax.set_ylabel("n. gains")
ax.set_xlabel("n. strains")

ax = axs[1]
sns.histplot(
    data=edf[~dupl_mask],
    y="n_losses",
    discrete=True,
    x="n_strains",
    ax=ax,
)
ax.plot([0, 111, 222], [0, 111, 0], "--", color="gray", alpha=0.3)
ax.set_ylabel("n. losses")
ax.set_xlabel("n. strains")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "nondupl_gainloss_vs_freq.png", dpi=300)
plt.show()

# %%
# single gain mask
sg_mask = edf["n_gains"] == 1
sg_mask &= edf["n_losses"] == 0
sg_mask &= edf["n_strains"] == 1
sg_mask &= edf["dupli"] == "no"
tgain_df = edf[sg_mask]

# single_loss mask
sl_mask = edf["n_gains"] == 0
sl_mask &= edf["n_losses"] == 1
sl_mask &= edf["n_strains"] == 222 - 1
sl_mask &= edf["dupli"] == "no"
tloss_df = edf[sl_mask]

# single gain loci
iso = "NZ_CP116085.1"
sdf = gdf.loc[tgain_df.index]
mask = sdf["locus"].str.startswith(iso)
gain_loci = sdf[mask].sort_values("locus")["locus"].str.removeprefix(iso + "_")

# single loss loci
other_iso = "NZ_CP098203.1"
sdf = gdf.loc[tloss_df.index]
mask = ~sdf["locus"].str.contains(iso)
loss_loci = sdf[mask]["locus"].str.split(other_iso + "_").str[1].str.split(" ").str[0]

# %%
# where do they start and end?


def parse_gbk(gbk_fname):
    with open(gbk_fname) as f:
        rec = SeqIO.read(f, "genbank")
    features = []
    for f in rec.features:
        if f.type == "CDS":
            features.append(
                {
                    "locus": f.qualifiers["locus_tag"][0],
                    "start": f.location.start,
                    "end": f.location.end,
                }
            )
    genome_len = len(rec.seq)
    return genome_len, pd.DataFrame(features).set_index("locus", verify_integrity=True)


fname = f"../../data/gbk/{iso}.gbk"
Li, gbk_df = parse_gbk(fname)
gain_starts = gbk_df.loc[gain_loci, "start"]

fname = f"../../data/gbk/{other_iso}.gbk"
Lo, gbk_df = parse_gbk(fname)
loss_starts = gbk_df.loc[loss_loci, "start"]

# %%


fig, axs = plt.subplots(2, 1, figsize=(8, 6))

ax = axs[0]
ax.hist(
    gain_starts,
    bins=np.arange(0, Li, 2000),
    cumulative=True,
    histtype="step",
)
ax.set_xlim(0, Li)
ax.set_xlabel(f"Genome {iso} position (bp)")
ax.set_ylabel("cumulative n. terminal gains")

ax = axs[1]
ax.hist(
    loss_starts,
    bins=np.arange(0, Lo, 2000),
    cumulative=True,
    histtype="step",
)
ax.set_xlim(0, Lo)
ax.set_xlabel(f"Genome {other_iso} position (bp)")
ax.set_ylabel("cumulative n. terminal losses")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "gainloss_pos.png", dpi=300)
plt.savefig(fig_fld / "gainloss_pos.pdf")
plt.show()

# %%
