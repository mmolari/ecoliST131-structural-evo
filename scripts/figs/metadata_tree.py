import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
from matplotlib.colors import to_hex


def parse_args():
    parser = argparse.ArgumentParser(description="Display metadata on core-genome tree")
    parser.add_argument("--tree", type=str)
    parser.add_argument("--metadata", type=str)
    parser.add_argument("--hist_fig", type=str)
    parser.add_argument("--tree_fig", type=str)
    return parser.parse_args()


def preprocess(df):
    # set host=wastewater when isolation_source=wastewater
    mask = df["isolation_source"] == "wastewater"
    df.loc[mask, "host"] = "wastewater"

    # date in bins
    df["year"] = pd.cut(
        df["year"],
        bins=[2002, 2008, 2012, 2016, 2020],
        include_lowest=True,
        labels=["2002-2008", "2009-2012", "2013-2016", "2017-2020"],
    )

    # order of location by frequency
    order_geo = df["geo_loc_name"].value_counts().index
    df["geo_loc_name"] = pd.Categorical(df["geo_loc_name"], order_geo, ordered=True)

    # order of host by frequency
    order_host = df["host"].value_counts().index
    df["host"] = pd.Categorical(df["host"], order_host, ordered=True)

    return df


def plot_histogram(df, svname):
    def addna(ax, series):
        n = series.isna().sum()
        ax.text(0.7, 0.1, f"N.A. = {n}", transform=ax.transAxes)

    fig, axs = plt.subplots(1, 3, figsize=(11, 4))

    ax = axs[0]
    sns.histplot(data=df, y="host", ax=ax)
    addna(ax, df["host"])
    ax.set_ylabel("host")

    ax = axs[1]
    sns.histplot(data=df, y="year", ax=ax)
    addna(ax, df["year"])
    ax.set_ylabel("year")

    ax = axs[2]
    sns.histplot(data=df, y="geo_loc_name", ax=ax)
    addna(ax, df["geo_loc_name"])
    ax.set_ylabel("location")

    sns.despine()
    plt.tight_layout()
    plt.savefig(svname)
    plt.close(fig)


def assign_colors(df):
    c_dict = {}

    cmap = plt.get_cmap("Dark2")
    cat = df["host"].dtype.categories
    c_dict["host"] = {v: to_hex(cmap(n)) for n, v in enumerate(cat)}
    c_dict["host"][np.nan] = "white"

    cmap = plt.get_cmap("viridis")
    cat = df["year"].dtype.categories
    c_dict["year"] = {v: to_hex(cmap(n / (len(cat) - 1))) for n, v in enumerate(cat)}
    c_dict["year"][np.nan] = "white"

    cmap = plt.get_cmap("tab20")
    cat = df["geo_loc_name"].dtype.categories
    c_dict["geo_loc_name"] = {
        v: to_hex(cmap((n * 2 % 20) + (n > 9))) for n, v in enumerate(cat)
    }
    c_dict["geo_loc_name"][np.nan] = "white"

    return c_dict


def plot_tree(df, tree_file, svname):
    # load tree
    tree = Phylo.read(tree_file, "newick")
    tree.ladderize()
    str_order = [x.name for x in tree.get_terminals()]
    N = len(str_order)

    # create color dictionary
    c_dict = assign_colors(df)

    # generate figure
    fig, axs = plt.subplots(1, 2, figsize=(8, 8), gridspec_kw={"width_ratios": [3, 1]})

    # draw tree
    ax = axs[0]
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: None)
    ax.grid(True, alpha=0.3)

    # add metadata markers
    xl, xr = ax.get_xlim()
    ddx = (xr - xl) * 0.05
    for i, k in enumerate(["host", "year", "geo_loc_name"]):
        D = df[k].to_dict()
        x = xr + ddx * np.ones(N) * (i - 2)
        y = np.arange(1, N + 1)
        c = [c_dict[k][D[s]] for s in str_order]
        ax.scatter(x, y, color=c, marker="o")
    ax.set_xlim(right=xr + ddx)
    ax.set_xlabel("branch length")
    ax.set_ylabel("n. isolates")
    for k in ["top", "right"]:
        ax.spines[k].set_visible(False)

    # legends
    ax = axs[1]
    ax.set_yticks([])
    E = []
    for i, k in enumerate(["host", "year", "geo_loc_name"]):
        elm = [
            mpl.lines.Line2D([], [], marker="o", ls="", color=c, label=l)
            for l, c in c_dict[k].items()
            if isinstance(l, str)
        ]
        E.append(elm)
    l0 = ax.legend(handles=E[0], title="host", loc="upper left")
    l1 = ax.legend(handles=E[1], title="collection year", loc="center left")
    l2 = ax.legend(handles=E[2], title="location", loc="lower left")
    ax.add_artist(l0)
    ax.add_artist(l1)

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    for k in ax.spines:
        ax.spines[k].set_visible(False)

    plt.tight_layout()
    plt.savefig(svname)
    plt.close(fig)


if __name__ == "__main__":
    args = parse_args()

    # load metadata
    df = pd.read_csv(args.metadata, index_col=0)

    # preprocess metadata
    df = preprocess(df)

    # plot histogram
    plot_histogram(df, args.hist_fig)

    # plot tree metadata
    plot_tree(df, args.tree, args.tree_fig)
