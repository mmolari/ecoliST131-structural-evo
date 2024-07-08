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
    parser.add_argument("--ann_df", type=str)
    parser.add_argument("--ann_gm", type=str)
    parser.add_argument("--ann_if", type=str)
    parser.add_argument("--ann_is", type=str)
    return parser.parse_args()


def add_ann_info(args, df):
    fnames = {
        "df": args.ann_df,
        "gm": args.ann_gm,
        "if": args.ann_if,
        "is": args.ann_is,
    }

    for k, fname in fnames.items():
        df2 = pd.read_csv(fname, index_col=0)
        df2 = df2["junction"].value_counts()
        df[f"{k}"] = df.index.map(df2)
        df[f"{k}"] = df[f"{k}"].fillna(0)
        df[f"{k}"] = df[f"{k}"].astype(int)

    return df


def load_df(args):
    df = pd.read_csv(args.junctions_stats, index_col=0)

    df2 = pd.read_csv(args.edge_pangenome, index_col=0)

    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    df["delta_L"] = df["max_length"] - df["min_length"]

    df = add_ann_info(args, df)

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
    np.random.seed(42)
    df["n_cat_wiggle"] = df["n_categories"] + np.random.uniform(-0.25, 0.25, len(df))
    df

    return df


def joint_distr(df, fig_fld):

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
    plt.close()


def fig_j_scatter_2(df, fig_fld):

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
    plt.close()


if __name__ == "__main__":

    args = parse_args()

    df = load_df(args)

    fig_fld = pathlib.Path(args.out_fld)
    fig_fld.mkdir(exist_ok=True, parents=True)

    joint_distr(df, fig_fld)
    fig_j_scatter(df, fig_fld)
    fig_j_scatter_2(df, fig_fld)
