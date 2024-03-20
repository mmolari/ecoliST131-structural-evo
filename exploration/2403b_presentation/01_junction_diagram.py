# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib


fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(parents=True, exist_ok=True)

# %%


def load_df():
    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
    df = pd.read_csv(fname, index_col=0)

    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_pangenome.csv"
    df2 = pd.read_csv(fname, index_col=0)

    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    df["delta_L"] = df["max_length"] - df["min_length"]
    df["acc_L"] = df["mean_length"] * 22
    return df


df = load_df()


def add_ann_info(df):
    fnames = {
        "df": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/defensefinder_real.csv",
        "gm": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/genomad_real.csv",
        "if": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/integronfinder_real.csv",
        "is": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/ISEScan_real.csv",
    }

    for k, fname in fnames.items():
        df2 = pd.read_csv(fname, index_col=0)
        df2 = df2["junction"].value_counts()
        df[f"{k}"] = df.index.map(df2)
        df[f"{k}"] = df[f"{k}"].fillna(0)
        df[f"{k}"] = df[f"{k}"].astype(int)

    return df


df = add_ann_info(df)
df = df[
    [
        # "n_iso",
        # "n_blocks",
        # "has_dupl",
        "n_categories",
        # "majority_category",
        # "singleton",
        # "cat_entropy",
        # "n_nodes",
        "min_length",
        "max_length",
        # "mean_length",
        "n_all_cores",
        # "core_left_length",
        # "core_right_length",
        "transitive",
        # "nonempty_acc_len",
        "nonempty_freq",
        "pangenome_len",
        # "pangenome_n_blocks",
        # "delta_L",
        # "acc_L",
        "df",
        "gm",
        "if",
        "is",
    ]
]
# remove transitive nodes
assert np.all(
    df["transitive"] == (df["n_categories"] == 1)
), "transitive means only one category"
assert np.all(df[df["transitive"]]["nonempty_freq"] == 0), "transitive means empty"
df = df[~df["transitive"]]
df
# %%
fig, axs = plt.subplots(
    2,
    2,
    figsize=(6, 5),
    sharex="col",
    sharey="row",
    gridspec_kw={"width_ratios": [4, 1], "height_ratios": [1, 4]},
)

ax = axs[1, 0]
sns.scatterplot(
    data=df,
    x="pangenome_len",
    y="n_categories",
    size=4,
    ax=ax,
    legend=False,
    alpha=0.4,
    marker="o",
    edgecolor="k",
    facecolor="none",
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("pangenome length (bp)")
ax.set_ylabel("n. categories")

kwargs = dict(
    stat="count",
    element="step",
    cumulative=False,
)

ax = axs[0, 0]
bins = np.histogram_bin_edges(np.log10(df["pangenome_len"]), bins=15)
sns.histplot(
    data=df,
    x="pangenome_len",
    color="gray",
    ax=ax,
    bins=bins,
    **kwargs,
)
ax.set_ylabel("n. junctions")

ax = axs[1, 1]
bins = np.histogram_bin_edges(np.log10(df["n_categories"]), bins=15)
sns.histplot(
    data=df,
    y="n_categories",
    color="gray",
    ax=ax,
    bins=bins,
    **kwargs,
)
ax.set_xlabel("n. junctions")

axs[0, 1].remove()

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"junct_diagram.png", dpi=150)
plt.show()

# %%
