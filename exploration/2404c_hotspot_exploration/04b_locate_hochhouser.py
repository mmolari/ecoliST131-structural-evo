# %%
import pandas as pd
import numpy as np
from Bio import Phylo
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt

# load tree
fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(fname, "newick")
strains = [l.name for l in tree.get_terminals()]

# load hochhauser hotspots
fname = "../../config/utils/hochhauser_SI.csv"
hhdf = pd.read_csv(fname, index_col=0)
cols = [
    "Gene symbol of upstream gene",
    "Gene symbol of downstream gene",
    "Number of genomes in which hotspot detected",
    "Number of genomes in which hotspot occupied",
    "Percentage of genomes in which hotspot occupied",
]
hhdf = hhdf[cols]

fname = "../../results/ST131_ABC/hotspots/asm20-100-5/hochhauser_junction_pos.csv"
pdf = pd.read_csv(fname, index_col=0)
pdf["gene"] = pdf.index.str.split("|").str[0]


fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
jdf = pd.read_csv(fname, index_col=0)
jdf
# %%
res = defaultdict(lambda: defaultdict(int))
for j, group in pdf.groupby("junction"):
    for iso in strains:
        subgroup = group[group["iso"] == iso]
        sgenes = subgroup["gene"].unique()
        for hs, row in hhdf.iterrows():
            gup = row["Gene symbol of upstream gene"]
            gdown = row["Gene symbol of downstream gene"]
            if gup in sgenes and gdown in sgenes:
                res[j][hs] += 1
res = pd.DataFrame(res).T
# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.set_aspect("equal", "box")
sns.heatmap(res, cmap="viridis")
# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
sns.histplot(data=pdf, y="gene", ax=ax)
plt.tight_layout()
plt.show()

# %%
fname = "../../results/ST131_ABC/hotspots/hochhauser_coregenes_position.csv"
df = pd.read_csv(fname, index_col=0)
# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
sns.countplot(data=df, y="type")
plt.xlim(200, 230)
plt.axvline(222, color="r", ls=":")
plt.tight_layout()
plt.show()

# %%
df.value_counts("type").loc["yjgR"]

# %%
found = set(df["type"].unique())

ok = 0
for hs, row in hhdf.iterrows():
    gup = row["Gene symbol of upstream gene"]
    gdown = row["Gene symbol of downstream gene"]
    if gup not in found:
        print(f"{hs} missing {gup}")
    if gdown not in found:
        print(f"{hs} missing {gdown}")
    if gup in found and gdown in found:
        print(f"{hs} oks!")
        ok += 1
print(ok)
# %%

set(jdf.index)
# %%
in_j = set(res.index) & set(jdf.index)
sdf = jdf.loc[list(in_j)]
# %%
sdf[sdf["n_categories"] >= 10]
# %%
