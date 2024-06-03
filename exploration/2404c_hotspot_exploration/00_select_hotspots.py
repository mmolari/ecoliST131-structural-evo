# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/f00")
fig_fld.mkdir(parents=True, exist_ok=True)

data_fld = pathlib.Path("data")
data_fld.mkdir(parents=True, exist_ok=True)


fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
df = pd.read_csv(fname, index_col=0)
fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_pangenome.csv"
df2 = pd.read_csv(fname, index_col=0)
df = pd.merge(df, df2, on="edge", validate="one_to_one")
mask = df["transitive"]
df = df[~mask]
df
# %%
min_len = 5e4
min_paths = 10

hs_mask = (df["pangenome_len"] >= min_len) & (df["n_categories"] >= min_paths)

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
sns.scatterplot(
    data=df,
    x="pangenome_len",
    y="n_categories",
    hue=hs_mask,
    legend=False,
    ax=ax,
    size=4,
)
plt.axvline(min_len, color="red", linestyle="--")
plt.axhline(min_paths, color="red", linestyle="--")
plt.xscale("log")
plt.yscale("log")
plt.text(0.1, 0.9, f"n. hotspots={hs_mask.sum()}", transform=plt.gca().transAxes)
plt.xlabel("pangenome length")
plt.ylabel("number of paths")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "hotspots.png")
plt.show()
# %%

# list of positive hotspots
hs = df[hs_mask].sort_values(
    ["pangenome_len", "n_categories"],
    ascending=False,
)
hs.to_csv(data_fld / "hotspots.csv")
hs
# %%
fname = "../../results/ST131_ABC/distances/summary-asm20-100-5.csv"
ddf = pd.read_csv(fname)
mask = ddf["si"] > ddf["sj"]
ddf = ddf[mask]
ddf
# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 3))
sns.histplot(data=ddf, x="core_div_filtered", bins=50, ax=ax)
ax.set_xlabel("filtered core genome divergence")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "core_div_filtered.png")
plt.show()

# %%
