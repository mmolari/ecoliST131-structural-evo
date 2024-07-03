# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib


fig_fld = pathlib.Path("figs/f02")
fig_fld.mkdir(parents=True, exist_ok=True)


def add_ann_info(df, dset):
    fnames = {
        "df": f"../../results/{dset}/annotations/junct_pos_asm20-100-5/defensefinder_real.csv",
        "gm": f"../../results/{dset}/annotations/junct_pos_asm20-100-5/genomad_real.csv",
        "if": f"../../results/{dset}/annotations/junct_pos_asm20-100-5/integronfinder_real.csv",
        "is": f"../../results/{dset}/annotations/junct_pos_asm20-100-5/ISEScan_real.csv",
    }

    for k, fname in fnames.items():
        df2 = pd.read_csv(fname, index_col=0)
        df2 = df2["junction"].value_counts()
        df[f"{k}"] = df.index.map(df2)
        df[f"{k}"] = df[f"{k}"].fillna(0)
        df[f"{k}"] = df[f"{k}"].astype(int)

    return df


def load_df(dset):
    fname = f"../../results/{dset}/backbone_joints/asm20-100-5/junctions_stats.csv"
    df = pd.read_csv(fname, index_col=0)

    fname = f"../../results/{dset}/backbone_joints/asm20-100-5/edge_pangenome.csv"
    df2 = pd.read_csv(fname, index_col=0)

    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    df["delta_L"] = df["max_length"] - df["min_length"]

    df = add_ann_info(df, dset)

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

    assert np.all(
        df["transitive"] == (df["n_categories"] == 1)
    ), "transitive means only one category"
    assert np.all(df[df["transitive"]]["nonempty_freq"] == 0), "transitive means empty"
    df = df[~df["transitive"]]
    df

    return df


# %%
for dset in ["ST131_ABC", "ST131_sub_BC", "ST131_sub_C", "ST131_sub_C2"]:
    df = load_df(dset)

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
                y="n_categories",
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
            y="n_categories",
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
            y="n_categories",
            alpha=0.8,
            color="firebrick",
            marker="o",
            size=8,
            legend=False,
            ax=ax,
            zorder=10,
        )
        ax.set_xlabel("pangenome length (bp)")
        ax.set_ylabel("n. distinct paths")
        ax.set_title(tt)

    plt.tight_layout()
    plt.savefig(fig_fld / f"MGE_{dset}.png", dpi=200)
    plt.show()


# %%

fig, axs = plt.subplots(
    2,
    2,
    figsize=(8, 8),
    sharex=True,
    sharey=True,
)
for nax, dset in enumerate(
    ["ST131_ABC", "ST131_sub_BC", "ST131_sub_C", "ST131_sub_C2"]
):
    df = load_df(dset)
    ax = axs.flatten()[nax]

    sns.scatterplot(
        data=df,
        x="pangenome_len",
        y="n_categories",
        color="k",
        alpha=0.2,
        size=4,
        ax=ax,
        legend=False,
    )

    ax.set_title(dset)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("pangenome length (bp)")
    ax.set_ylabel("n. distinct paths")
    ax.set_xlim(10, 1e6)
    ax.set_ylim(1, 222)
    ax.grid(axis="both", which="major", linewidth=0.5, alpha=0.2)
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"junction_survey.png", dpi=200)
plt.show()

# %%

dset_color = {
    "ST131_ABC": "C0",
    "ST131_sub_BC": "C1",
    "ST131_sub_C": "C2",
    "ST131_sub_C2": "C3",
}

fig, axs = plt.subplots(1, 2, figsize=(8, 4), sharey=True)

bins_cat = np.logspace(0, 3, 200)
bins_len = np.logspace(1, 6, 200)

for dset in ["ST131_ABC", "ST131_sub_BC", "ST131_sub_C", "ST131_sub_C2"]:
    df = load_df(dset)

    ax = axs[0]
    sns.histplot(
        data=df,
        x="n_categories",
        color=dset_color[dset],
        ax=ax,
        element="step",
        fill=False,
        bins=bins_cat,
        stat="count",
        cumulative=True,
        label=dset,
    )

    ax = axs[1]
    sns.histplot(
        data=df,
        x="pangenome_len",
        color=dset_color[dset],
        ax=ax,
        element="step",
        fill=False,
        bins=bins_len,
        stat="count",
        label=dset,
        cumulative=True,
    )
axs[0].set_xscale("log")
axs[1].set_xscale("log")
axs[0].set_xlabel("n. distinct paths")
axs[1].set_xlabel("pangenome length (bp)")
axs[0].set_ylabel("n. junctions")
plt.legend()
plt.tight_layout()
plt.savefig(fig_fld / "junct_cumdistr.png", dpi=200)
plt.show()

# %%
