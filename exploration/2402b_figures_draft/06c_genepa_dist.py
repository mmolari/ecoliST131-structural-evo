# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
import pathlib
from collections import Counter
from Bio import AlignIO
import numpy as np

fig_fld = pathlib.Path("figs/f06")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/distances/summary-asm20-100-5.csv"
df = pd.read_csv(fname)
mask = df["si"] > df["sj"]
df = df[mask]


with open(
    "../../results/ST131_ABC/pangraph/asm20-100-5-alignment/filtered_corealignment_info_size.json"
) as f:
    aln_info = json.load(f)
core_aln_size = aln_info["core aln size"]
red_aln_size = aln_info["polished aln size"]

df["filt_aln_SNPs"] = df["core_div_filtered"] * red_aln_size

cols = [
    "si",
    "sj",
    "core_div_naive",
    "mash_dist",
    "private seq. (bp)",
    "shared seq. (bp)",
    "n. breakpoints",
    "part. entropy",
    "n. blocks",
    "core_div_filtered",
    "edge_PA",
    "edge_PA_reduced",
    "edge_sharing",
    "block_PA",
    "block_sharing",
    "acc_block_PA",
]

strains = df["si"].unique()


# load gene clusters info
gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"

with open(gene_cluster_json, "r") as f:
    gdf = pd.DataFrame(json.load(f))
gdf["divers"] = gdf["divers"].astype(float)
gdf["count"] = gdf["count"].astype(int)
gdf["geneLen"] = gdf["geneLen"].astype(int)
gdf["event"] = gdf["event"].astype(int)


# get gene count
pth = pathlib.Path("../../results/ST131_ABC/panx/gc_loci")
gene_counts = {}
for gid in gdf["geneId"]:
    gc = pd.read_csv(pth / f"genecl_{gid}.csv")
    ct = Counter(gc["iso"])
    gene_counts[gid] = dict(ct)
gene_counts = pd.DataFrame(gene_counts)
gene_counts = gene_counts.loc[strains]
gene_counts.fillna(0, inplace=True)
gene_counts = gene_counts.astype(int)
gene_counts = (gene_counts > 0).astype(int)

sns.heatmap(gene_counts)
plt.show()


# %%

mask = df["si"] == "NC_022648.1"
mask |= df["sj"] == "NC_022648.1"
sdf = df[~mask].copy()


def grr(x, y):
    return np.dot(x, y) / max(x.sum(), y.sum())


sdf["GRR"] = sdf.apply(
    lambda row: grr(
        gene_counts.loc[row["si"]].to_numpy(), gene_counts.loc[row["sj"]].to_numpy()
    ),
    axis=1,
)

# %%

fig, ax = plt.subplots(figsize=(6, 4))
sns.histplot(data=sdf, x="core_div_filtered", y="GRR", bins=50)

# plot median for each bin
bins = np.linspace(0, sdf["core_div_filtered"].max(), 20)
sdf["bin"] = pd.cut(sdf["core_div_filtered"], bins=bins)
medians = sdf.groupby("bin")["GRR"].median()
plt.plot(np.diff(bins) / 2 + bins[:-1], medians, "k:")


plt.xlabel("core-genome divergence (filtered)")
plt.xticks(np.arange(0, 1.7e-4, 0.5e-4))
plt.ylabel("GRR")
plt.tight_layout()
sns.despine()
plt.savefig(fig_fld / "grr_vs_divergence.png")
plt.savefig(fig_fld / "grr_vs_divergence.pdf")
plt.show()
# %%
