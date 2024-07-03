# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import numpy as np
import json

import pypangraph as pp
import utils as ut

from collections import defaultdict
from Bio import Phylo


dset = "ST131_ABC"
fig_fld = pathlib.Path(f"figs/n00a/{dset}")
fig_fld.mkdir(exist_ok=True, parents=True)


def svfig(name):
    # plt.savefig(fig_fld / f"{name}.svg")
    plt.savefig(fig_fld / f"{name}.png", facecolor="white", dpi=150)


fld = pathlib.Path(f"../../results/{dset}")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_edge_len = fld / "backbone_joints/asm20-100-5/edge_len.csv"
pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

with open(pos_file) as f:
    jp = json.load(f)
df = pd.read_csv(df_file, index_col=0)
df_el = pd.read_csv(df_edge_len, index_col=0)
df["delta_len"] = df["max_length"] - df["min_length"]


def extend_df(df, df_el):
    new_df = {}
    # FG = full graph
    # frequency of edges
    new_df["edge_freq_FG"] = df_el.notna().sum(axis=0) / df_el.shape[0]
    new_df["nonempty_freq_FG"] = (df_el > 0).sum(axis=0) / df_el.notna().sum(axis=0)
    new_df["n_iso_FG"] = df_el.notna().sum(axis=0)
    new_df["nonempty_acc_len_FG"] = df_el.sum(axis=0) / (df_el > 0).sum(axis=0)
    new_df = pd.DataFrame(new_df)
    new_df.index.name = "edge"
    df = df.join(new_df, how="outer").sort_values("n_iso_FG", ascending=False)
    return df


df = extend_df(df, df_el)


def edge_cat(freq):
    if freq < 0.1:
        return "rare"
    elif freq < 0.9:
        return "intermediate"
    elif freq < 1.0:
        return "common"
    elif freq == 1.0:
        return "backbone"
    else:
        raise ValueError(f"invalid freq: {freq}")


df["edge_cat"] = df["edge_freq_FG"].apply(edge_cat)

df

# %% stats

mask = df["nonempty_acc_len_FG"].isna()
# mask = df["edge_freq_FG"].isna()
df[mask]
# %% transitive

mask = df["n_categories"] == 1
df[mask].to_csv(f"data/{dset}/transitive.csv")
df[mask]

# %% edge frequency

freqs = df["n_iso_FG"].value_counts()

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
sns.barplot(x=freqs.index, y=freqs.values, ax=ax, color="gray", log=True)
ax.bar_label(ax.containers[0], fmt="%d")
ax.set_title("edge frequency")
ax.set_xlabel("n. isolates")
ax.set_ylabel("n. edges")
plt.tight_layout()
svfig("edge_freq_FG")
plt.show()

# %% edge length and occupation frequency, stratified by category

g = sns.displot(
    data=df,
    x="nonempty_freq_FG",
    y="nonempty_acc_len_FG",
    kind="hist",
    bins=50,
    col="edge_cat",
    height=4,
    aspect=0.8,
    log_scale=(False, True),
)
for ax in g.axes.flat:
    ax.set_xlabel("edge occupation frequency")
    ax.set_ylabel("occupied edge average length (bp)")
    title = ax.get_title().split("=")[1].strip()
    n_points = df["edge_cat"] == title
    n_points &= df["nonempty_acc_len_FG"].notna()
    title += f" (n={n_points.sum()})"
    ax.set_title(title)
    ax.grid(alpha=0.3)

plt.tight_layout()
svfig("edge_cats_FG")
plt.show()
# %%

mask = df["n_iso"] > 50
mask &= df["n_categories"] > 1
sdf = df[mask].copy()

g = sns.JointGrid(
    data=sdf, x="nonempty_freq", y="nonempty_acc_len", marginal_ticks=True
)
g.ax_joint.set(yscale="log")
g.plot_joint(sns.scatterplot, alpha=0.1)
g.plot_marginals(sns.histplot, element="step", bins=50)
g.ax_joint.set_xlabel("edge occupation frequency")
g.ax_joint.set_ylabel("occupied edge average length (bp)")
g.ax_joint.grid(alpha=0.3)
# plt.tight_layout()
svfig("edge_freq_len_joint")
plt.show()
# %%

f1, f2 = "nonempty_freq", "nonempty_freq_FG"
l1, l2 = "nonempty_acc_len", "nonempty_acc_len_FG"

sdf["weird"] = np.abs(np.log10(sdf[l2] / sdf[l1])) > 0.1

for svname, kwargs in [["weird", {"hue": "weird"}], ["comparison", {}]]:
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    ax = axs[0]
    sns.scatterplot(data=sdf, x=f1, y=f2, ax=ax, alpha=0.1, **kwargs)
    ax.grid(alpha=0.3)
    ax.set_xlabel("edge occupation frequency (refined)")
    ax.set_ylabel("edge occupation frequency (full graph)")
    ax = axs[1]
    sns.scatterplot(data=sdf, x=l1, y=l2, ax=ax, alpha=0.1, **kwargs)
    ax.grid(alpha=0.3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("occupied edge average length (bp) (full graph)")
    ax.set_xlabel("occupied edge average length (bp) (refined)")
    plt.tight_layout()
    svfig(svname)
    plt.show()


sdf[sdf["weird"]].to_csv(f"data/{dset}/weird.csv")
sdf[sdf["weird"]]


# %%

g = sns.jointplot(
    data=sdf,
    x="delta_len",
    y="n_categories",
    kind="scatter",
    marginal_kws={"log_scale": True},
    joint_kws={"alpha": 0.2},
    space=0,
)
g.ax_joint.set_xlabel("len max - len min (bp)")
g.ax_joint.set_ylabel("n. distinct paths")
g.ax_joint.grid(alpha=0.3)
g.ax_joint.set_ylim(bottom=1)

# plt.tight_layout()
svfig("delta_len_vs_ncats")
plt.show()

# %%
