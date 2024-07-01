import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_fld", type=str, default="figs")
    parser.add_argument("--junctions_stats", type=str)
    parser.add_argument("--edge_pangenome", type=str)
    return parser.parse_args()


def load_df(args):
    df = pd.read_csv(args.junctions_stats, index_col=0)

    df2 = pd.read_csv(args.edge_pangenome, index_col=0)

    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    df["delta_L"] = df["max_length"] - df["min_length"]

    df = df[
        [
            "n_iso",
            "n_categories",
            "min_length",
            "max_length",
            "n_all_cores",
            "transitive",
            "nonempty_freq",
            "pangenome_len",
        ]
    ]

    assert np.all(
        df["transitive"] == (df["n_categories"] == 1)
    ), "transitive means only one category"
    assert np.all(df[df["transitive"]]["nonempty_freq"] == 0), "transitive means empty"
    df = df[~df["transitive"]]
    np.random.seed(42)
    df["n_cat_wiggle"] = df["n_categories"] + np.random.uniform(-0.25, 0.25, len(df))
    df

    return df


def fig_occ_freq(df, fig_fld):

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
    ax.set_xlabel("pangenome length (bp)")
    ax.set_ylabel("n. categories")
    plt.tight_layout()
    plt.savefig(fig_fld / "diag_occ_freq.pdf")
    plt.savefig(fig_fld / "diag_occ_freq.svg")
    plt.close()


def fig_j_scatter(df, fig_fld):

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
        facecolor="none",
        edgecolor=color,
        size=4,
        alpha=0.8,
        legend=False,
        rasterized=True,
        ax=ax,
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
    plt.close()


def fig_n_paths(df, fig_fld):

    sdf = df["n_categories"] - 1
    sdf = pd.DataFrame(sdf)

    Nmax = sdf["n_categories"].max()
    discrete_bins_cat = np.arange(Nmax + 2) - 0.5

    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax = axs[0]
    sns.histplot(
        sdf,
        cumulative=True,
        element="step",
        fill=False,
        bins=discrete_bins_cat,
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
        bins=discrete_bins_cat,
        binrange=(0, 150),
        ax=ax,
    )

    ax.set_ylabel("cumulative sum")
    ax.set_xlabel("n. path categories - 1")
    ax.set_xlim(right=150)
    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / "cumulative_categories.pdf")
    plt.savefig(fig_fld / "cumulative_categories.svg")
    plt.close()


if __name__ == "__main__":

    args = parse_args()

    df = load_df(args)

    fig_fld = pathlib.Path(args.out_fld)
    fig_fld.mkdir(exist_ok=True, parents=True)

    fig_occ_freq(df, fig_fld)
    fig_j_scatter(df, fig_fld)
    fig_n_paths(df, fig_fld)
