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
fig_fld = pathlib.Path(f"figs/n0/{dset}")
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
df

# %% stats

mask = df["nonempty_acc_len_FG"].isna()
# mask = df["edge_freq_FG"].isna()
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
    # ax.set_yscale("symlog", linthresh=1)
    # ax.set_ylim(bottom=0)
    # ax.set_yscale("log")
    # ax.set_xlim(left=0)
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

sdf = df[df["n_iso"] > 50]

g = sns.JointGrid(
    data=sdf, x="nonempty_freq", y="nonempty_acc_len", marginal_ticks=True
)
g.ax_joint.set(yscale="log")
g.plot_joint(sns.scatterplot, alpha=0.1, color="#03012d")
g.plot_marginals(sns.histplot, element="step", color="#03012d", bins=50)
g.ax_joint.set_xlabel("edge occupation frequency")
g.ax_joint.set_ylabel("occupied edge average length (bp)")
g.ax_joint.grid(alpha=0.3)
# set title
title = f"n. isolates = {sdf.shape[0]}"
plt.tight_layout()
svfig("edge_freq_len_joint")
plt.show()
# %%
freq_1 = (df_el > 0).sum(axis=0) / df_el.notna().sum(axis=0)
len_1 = df_el.sum(axis=0) / (df_el > 0).sum(axis=0)
freq_2 = df["nonempty_freq"]
len_2 = df["nonempty_acc_len"]

sdf = pd.DataFrame(
    {
        "freq_1": freq_1,
        "len_1": len_1,
        "freq_2": freq_2,
        "len_2": len_2,
        "transitive": df["transitive"],
    },
    index=freq_1.index,
)

# only_iso = df[df["n_iso"] > 50].index
# sdf = sdf.loc[only_iso]

sdf["weird"] = sdf["freq_1"] - sdf["freq_2"] > 0.5
# sdf["weird"] = sdf["len_2"] / sdf["len_1"] > 200

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
sns.scatterplot(data=sdf, x="freq_1", y="freq_2", ax=axs[0], alpha=0.1, hue="weird")
ax = axs[1]
sns.scatterplot(data=sdf, x="len_1", y="len_2", ax=ax, alpha=0.1, hue="weird")
ax.grid(alpha=0.3)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax
plt.tight_layout()
plt.show()


# %%
ssdf = df[sdf["weird"]]
ssdf
# %%
sdf[sdf["len_2"].isna()]["len_1"].hist(bins=np.logspace(0, 5, 50))
plt.xscale("log")
plt.show()
# %%
