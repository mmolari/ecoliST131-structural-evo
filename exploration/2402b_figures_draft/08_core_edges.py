# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/f08")
fig_fld.mkdir(parents=True, exist_ok=True)

fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_count.csv"
dfe = pd.read_csv(fname, index_col=0)


def load_df():
    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
    df = pd.read_csv(fname, index_col=0)

    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_pangenome.csv"
    df2 = pd.read_csv(fname, index_col=0)

    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    df["delta_L"] = df["max_length"] - df["min_length"]
    return df


df = load_df()
df
# %%
dfe["transitive"] = df["transitive"].astype(bool)
mask = dfe["transitive"].isna()
dfe.loc[mask, "transitive"] = True
dfe["transitive"] = dfe["transitive"].astype(bool)
# %%

fig, ax = plt.subplots(1, 1, figsize=(7, 4))

Ns = np.sort(dfe["count"].unique())
M = Ns.max()
nl = Ns[Ns < (M / 2)].max()
nh = Ns[Ns >= (M / 2)].min()
Ns = list(range(1, nl + 2)) + list(range(nh, M + 1))

for i, N in enumerate(Ns):
    mask = dfe["count"] == N
    smask = mask & (~dfe["transitive"])
    k = smask.sum()
    ax.bar(i - 0.2, k, 0.4, color="dimgray")
    if k > 0:
        ax.text(i - 0.2, k + 0.4, k, ha="center", va="bottom", color="dimgray")
    smask = mask & dfe["transitive"]
    k = smask.sum()
    ax.bar(i + 0.2, k, 0.4, color="silver")
    if k > 0:
        ax.text(i + 0.2, k + 0.4, k, ha="center", va="bottom", color="silver")

#  make legend
ax.bar(0, 0, color="dimgray", label="non-transitive")
ax.bar(0, 0, color="silver", label="transitive")
ax.legend()

ax.set_xticks(range(len(Ns)))
ax.set_xticklabels(Ns)
ax.set_xlabel("Number of isolates with the edge")
ax.set_ylabel("Number of edges")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "edge_count.svg")
plt.savefig(fig_fld / "edge_count.png")
plt.show()
# %%
