# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib


fig_fld = pathlib.Path("figs/f03")
fig_fld.mkdir(parents=True, exist_ok=True)

# %%


def load_df():
    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
    df = pd.read_csv(fname, index_col=0)

    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_pangenome.csv"
    df2 = pd.read_csv(fname, index_col=0)

    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    df["delta_L"] = df["max_length"] - df["min_length"]
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
        "n_iso",
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
np.random.seed(42)
df["n_cat_wiggle"] = df["n_categories"] + np.random.uniform(-0.25, 0.25, len(df))
df
# %%

fig, ax = plt.subplots(1, 1, figsize=(7, 4))
norm = mpl.colors.Normalize(vmin=0, vmax=1)
mapp = mpl.cm.ScalarMappable(norm=norm, cmap="coolwarm")
sns.scatterplot(
    data=df,
    x="pangenome_len",
    y="n_cat_wiggle",
    alpha=0.3,
    hue="nonempty_freq",
    hue_norm=norm,
    palette="coolwarm",
    size=4,
    ax=ax,
    legend=False,
)
plt.colorbar(mapp, ax=ax, label="occupation frequency")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("local pangenome length (bp)")
ax.set_ylabel("n. distinct paths")
plt.tight_layout()
plt.savefig(fig_fld / "diag_occ_freq.pdf")
plt.savefig(fig_fld / "diag_occ_freq.svg")
plt.show()

# %%


fig, axs = plt.subplots(2, 2, figsize=(10, 7.2), sharex=True, sharey=True)
for lab, ax, bw, col, cm, tt in [
    ("gm", axs[0, 0], 0.5, "C0", "Blues", "prophages"),
    ("if", axs[0, 1], 0.5, "C1", "Oranges", "integrons"),
    ("df", axs[1, 0], 0.5, "C2", "Greens", "defense systems"),
    ("is", axs[1, 1], 0.5, "C3", "Purples", "insertion sequences"),
]:
    mask = df[lab] > 0

    if bw > 0:
        sns.kdeplot(
            data=df[mask],
            x="pangenome_len",
            y="n_cat_wiggle",
            fill=True,
            levels=5,
            # linewidths=3,
            cmap=cm,
            log_scale=(True, True),
            ax=ax,
            bw_adjust=bw,
        )
    sns.scatterplot(
        data=df[~mask],
        x="pangenome_len",
        y="n_cat_wiggle",
        alpha=0.5,
        edgecolors="k",
        facecolors="none",
        size=4,
        ax=ax,
        legend=False,
        zorder=9,
    )

    sns.scatterplot(
        data=df[mask],
        x="pangenome_len",
        y="n_cat_wiggle",
        alpha=0.8,
        color="firebrick",
        marker="o",
        size=8,
        legend=False,
        ax=ax,
        zorder=10,
    )
    ax.set_xlabel("local pangenome length (bp)")
    ax.set_ylabel("n. distinct paths")
    ax.set_title(tt)

plt.tight_layout()
plt.savefig(fig_fld / f"diag_suppl.pdf")
plt.savefig(fig_fld / f"diag_suppl.svg")
plt.show()

# %%


fig, axs = plt.subplots(2, 2, figsize=(10, 7.2), sharex=True, sharey=True)

ax = axs[0, 0]
norm = mpl.colors.Normalize(vmin=0, vmax=1)
mapp = mpl.cm.ScalarMappable(norm=norm, cmap="coolwarm")
sns.scatterplot(
    data=df,
    x="pangenome_len",
    y="n_cat_wiggle",
    alpha=0.4,
    hue="nonempty_freq",
    hue_norm=norm,
    palette="coolwarm",
    size=4,
    ax=ax,
    legend=False,
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("local pangenome length (bp)")
ax.set_ylabel("n. distinct paths")


for lab, ax, bw, col, cm, tt in [
    ("is", axs[0, 1], 0.5, "C0", "Blues", "insertion sequences"),
    ("df", axs[1, 0], 0.5, "C2", "Greens", "defense systems"),
    ("gm", axs[1, 1], 0.5, "C4", "Purples", "prophages"),
]:
    mask = df[lab] > 0

    sns.kdeplot(
        data=df[mask],
        x="pangenome_len",
        y="n_cat_wiggle",
        fill=True,
        levels=8,
        cmap=cm,
        log_scale=(True, True),
        ax=ax,
        bw_adjust=0.5,
    )
    sns.scatterplot(
        data=df[~mask],
        x="pangenome_len",
        y="n_cat_wiggle",
        alpha=0.5,
        edgecolors="k",
        facecolors="none",
        size=4,
        ax=ax,
        legend=False,
        zorder=9,
    )

    sns.scatterplot(
        data=df[mask],
        x="pangenome_len",
        y="n_cat_wiggle",
        alpha=0.8,
        color="firebrick",
        marker="o",
        size=4,
        legend=False,
        ax=ax,
        zorder=10,
    )
    ax.set_xlabel("local pangenome length (bp)")
    ax.set_ylabel("n. distinct paths")
    ax.set_title(tt, color=col, fontsize=14)

plt.tight_layout()

cax = fig.add_axes([0.08, 0.74, 0.01, 0.15])
plt.colorbar(mapp, cax=cax, label="occupation frequency")

plt.savefig(fig_fld / f"diag_main.pdf")
plt.savefig(fig_fld / f"diag_main.svg")
plt.show()

# %%

for lab, col, tt in [
    ("is", "C0", "insertion sequences"),
    ("if", "C3", "integrons"),
    ("df", "C2", "defense systems"),
    ("gm", "C4", "prophages"),
]:
    mask = df[lab] > 0
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
        y="n_cat_wiggle",
        color="lightgrey",
        size=4,
        ax=ax,
        legend=False,
        rasterized=True,
        # label="total",
    )
    sns.scatterplot(
        data=df[mask],
        x="pangenome_len",
        y="n_cat_wiggle",
        color=col,
        size=4,
        ax=ax,
        legend=False,
        label=tt,
    )
    ax.legend(loc="upper left")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("local pangenome length (bp)")
    ax.set_ylabel("n. distinct paths")

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
        color="lightgray",
        ax=ax,
        bins=bins,
        **kwargs,
    )
    sns.histplot(
        data=df[mask],
        x="pangenome_len",
        color=col,
        ax=ax,
        bins=bins,
        **kwargs,
    )
    ax.set_ylabel("n. junctions")

    ax = axs[1, 1]
    bins = np.histogram_bin_edges(np.log10(df["n_categories"]), bins=15)
    sns.histplot(
        data=df,
        y="n_cat_wiggle",
        color="lightgray",
        ax=ax,
        bins=bins,
        **kwargs,
    )
    sns.histplot(
        data=df[mask],
        y="n_cat_wiggle",
        color=col,
        ax=ax,
        bins=bins,
        **kwargs,
    )
    ax.set_xlabel("n. junctions")

    axs[0, 1].remove()

    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / f"diag_suppl_{lab}.pdf", dpi=200)
    plt.savefig(fig_fld / f"diag_suppl_{lab}.svg", dpi=200)
    plt.show()

# %%

color = "grey"

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
    y="n_cat_wiggle",
    # color=color,
    # facecolor=color,
    facecolor="none",
    edgecolor=color,
    size=4,
    alpha=0.8,
    legend=False,
    rasterized=True,
    ax=ax,
)
# sns.kdeplot(
#     data=df,
#     x="pangenome_len",
#     y="n_cat_wiggle",
#     fill=False,
#     levels=5,
#     cmap="Greys",
#     vmin=-0.5,
#     log_scale=(True, True),
#     zorder=-5,
#     ax=ax,
# )
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("local pangenome length (bp)")
ax.set_ylabel("n. distinct paths")

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
    color=color,
    ax=ax,
    bins=bins,
    **kwargs,
)
ax.set_ylabel("n. junctions")

ax = axs[1, 1]
bins = np.histogram_bin_edges(np.log10(df["n_cat_wiggle"]), bins=15)
sns.histplot(
    data=df,
    y="n_cat_wiggle",
    color=color,
    ax=ax,
    bins=bins,
    **kwargs,
)
ax.set_xlabel("n. junctions")

axs[0, 1].remove()

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"simple_scatter.pdf", dpi=300)
plt.savefig(fig_fld / f"simple_scatter.svg", dpi=300)
plt.show()

# %%

fig, axs = plt.subplots(3, 1, figsize=(4, 8), sharey=True, sharex=True)

for ax, lab, col, tt in [
    (axs[0], "is", "C0", "insertion sequences"),
    (axs[1], "gm", "C4", "prophages"),
    (axs[2], "df", "C2", "defense systems"),
]:
    mask = df[lab] > 0
    sns.scatterplot(
        data=df,
        x="pangenome_len",
        y="n_cat_wiggle",
        edgecolor="gray",
        facecolor="none",
        size=4,
        ax=ax,
        legend=False,
        rasterized=True,
    )
    sns.scatterplot(
        data=df[mask],
        x="pangenome_len",
        y="n_cat_wiggle",
        facecolor=col,
        edgecolor=col,
        alpha=0.5,
        size=4,
        ax=ax,
        legend=False,
        rasterized=True,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("local pangenome length (bp)")
    ax.set_ylabel("n. distinct paths")
    ax.set_title(tt)

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"simple_scatter2.pdf", dpi=300)
plt.savefig(fig_fld / f"simple_scatter2.svg", dpi=300)
plt.show()
# %%

sdf = df["n_categories"] - 1
sdf = pd.DataFrame(sdf)
fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
ax = axs[0]
sns.histplot(
    sdf,
    cumulative=True,
    element="step",
    fill=False,
    discrete=True,
    ax=ax,
    legend=False,
    binrange=(0, 150),
)
ax.set_ylabel("cumulative distribution")

ax = axs[1]
sns.histplot(
    x=sdf["n_categories"],
    cumulative=True,
    weights=sdf["n_categories"],
    element="step",
    fill=False,
    discrete=True,
    binrange=(0, 150),
    ax=ax,
)

ax.set_ylabel("cumulative sum")
ax.set_xlabel("n. distinct paths - 1")
ax.set_xlim(right=150)
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "cumulative_categories.pdf")
plt.savefig(fig_fld / "cumulative_categories.svg")
plt.show()
# %%

cms = sdf.sort_values(by="n_categories").cumsum()
mask = cms > cms.max() * 0.5
sdf[mask].min()
# %%
