import numpy as np
import pandas as pd
import pathlib as pthl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from Bio import Phylo
import argparse

# assign a color to each continent
continent_color = {
    "Africa": "tab:purple",
    "Asia": "tab:red",
    "Europe": "tab:blue",
    "North America": "tab:green",
    "Oceania": "tab:orange",
    "South America": "tab:brown",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Plot metadata")
    parser.add_argument("--metadata_csv", type=str)
    parser.add_argument("--alleles_csv", type=str)
    parser.add_argument("--coregenome_tree", type=str)
    parser.add_argument("--outdir", type=str)
    return parser.parse_args()


def country_barplot(df, svname, **kwargs):
    continent_count = df["continent"].value_counts()

    # make legend
    legend_elements = []
    for k, v in continent_color.items():
        legend_elements.append(
            mpl.lines.Line2D(
                [0],
                [0],
                color=v,
                label=k + f" (n={continent_count[k]})",
                marker="s",
                ls="",
            )
        )

    # count entries per country
    counts = df[["geo_loc_name", "continent"]].value_counts()
    N = len(counts)

    fig, ax = plt.subplots(1, 1, figsize=(2 + N * 0.25, 4.5))
    # draw colored barplot
    for idx, row in counts.iteritems():
        ctr, cont = idx
        ax.bar(ctr, row, color=continent_color[cont])
    # set labels
    ax.set_xlabel("Country")
    ax.set_ylabel("Number of isolates")
    # rotate xticks
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center")
    ax.legend(handles=legend_elements, loc="upper right", title="Continent")
    sns.despine()
    plt.tight_layout()
    plt.savefig(svname, **kwargs)
    plt.close(fig)


def collection_date_hist(df, svname, **kwargs):
    yrs = df["year"]
    M, m = yrs.max(), yrs.min()
    bins = np.arange(m, M + 3) - 0.5
    # licks every 5 years for round years
    fig, ax = plt.subplots(
        1,
        1,
    )
    ax.hist(yrs, bins=bins, edgecolor="black", color="lightgrey")
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of isolates")
    xticks = [y for y in range(int(m), int(M) + 2) if y % 5 == 0]
    ax.set_xticks(xticks)
    sns.despine()
    plt.tight_layout()
    plt.savefig(svname, **kwargs)
    plt.close(fig)


def __mark(ax, x, y, c, **kwargs):
    ax.scatter(x, y, marker="s", edgecolor="lightgray", color=c, **kwargs)


def __create_legend(color_dictionaries):
    L = len(color_dictionaries)
    fig, axs = plt.subplots(1, L, figsize=(L * 1.5, 5))
    for i, (name, color_dict) in enumerate(color_dictionaries.items()):
        ax = axs[i]
        for k, v in color_dict.items():
            __mark(ax, [], [], c=v, label=k)
        ax.legend(
            loc="upper left",
            title=name,
            frameon=False,
        )
        ax.axis("off")
    plt.tight_layout()
    plt.subplots_adjust(wspace=0)
    return fig, axs


def metadata_tree(df_met, df_all, tree, svname_plot, svname_leg, **kwargs):
    alleles = list(df_all.columns)
    df = pd.concat([df_met, df_all], axis=1)

    N = len(tree.get_terminals())
    fig, axs = plt.subplots(
        1,
        2,
        figsize=(5, 1 + 0.1 * N),
        sharey=True,
        gridspec_kw={"width_ratios": [1, 0.3]},
    )

    # plot tree
    ax = axs[0]
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center")

    # plot metadata
    ax = axs[1]
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # assign colors to the 12 most common countries
    countries = df["geo_loc_name"].value_counts()
    top_countries = countries.index[:10]
    cmap = plt.cm.get_cmap("tab10")
    country_color = {k: cmap(i) for i, k in enumerate(top_countries)}

    # assign color to year
    ym, yM = int(df["year"].min()), int(df["year"].max())
    cmap = plt.cm.get_cmap("Blues")
    yr_range = range(ym, yM + 1)
    year_color = {y: cmap(i / (len(yr_range) - 1)) for i, y in enumerate(yr_range)}

    # assign color to isolation
    isolations = df["isolation_source"].value_counts()
    cmap = plt.cm.get_cmap("tab10")
    isolation_color = {i: cmap(n) for n, i in enumerate(isolations.index)}

    # assign color to alleles
    allele_colors = {}
    for al in alleles:
        all_ct = df[al].value_counts()
        cmap = plt.cm.get_cmap("tab10")
        allele_colors[al] = {int(i): cmap(n) for n, i in enumerate(all_ct.index)}

    xticks = {"country": 0, "source": 1, "year": 2}
    xticks |= {al: n + 4 for n, al in enumerate(alleles)}
    for i, name in enumerate([l.name for l in tree.get_terminals()]):
        y = i + 1
        # get metadata
        row = df.loc[name]
        country = row["geo_loc_name"]
        year = row["year"]
        isol = row["isolation_source"]

        if pd.notna(country):
            if country in top_countries:
                c = country_color[country]
            else:
                c = "white"
            x = xticks["country"]
            __mark(ax, x, y, c)
        if pd.notna(isol):
            x = xticks["source"]
            __mark(ax, x, y, c=isolation_color[isol])
        if pd.notna(year):
            x = xticks["year"]
            __mark(ax, x, y, c=year_color[int(year)])
        for al in alleles:
            al_r = row[al]
            if pd.notna(al_r):
                x = xticks[al]
                __mark(ax, x, y, c=allele_colors[al][int(al_r)])

    ax.set_xticks(list(xticks.values()))
    ax.set_xticklabels(list(xticks.keys()), rotation=90)

    for ax in axs:
        for i in range(N):
            ax.axhline(i + 1, color="grey", ls="--", lw=0.5, zorder=-1)

    # add color legents

    plt.tight_layout()
    plt.subplots_adjust(wspace=0)
    plt.savefig(svname_plot, **kwargs)
    plt.close(fig)

    color_dicts = {
        "country": country_color,
        "source": isolation_color,
        "year": year_color,
    } | allele_colors

    fig, axs = __create_legend(color_dicts)
    plt.savefig(svname_leg, **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    # load input data
    args = parse_args()
    df_met = pd.read_csv(args.metadata_csv, index_col=0)
    df_all = pd.read_csv(args.alleles_csv, index_col=0)
    tree = Phylo.read(args.coregenome_tree, "newick")
    svfld = pthl.Path(args.outdir)
    svfld.mkdir(exist_ok=True, parents=True)

    # barplot by country
    country_barplot(df_met, svfld / "country_barplot.png", dpi=300)

    # collection date histogram
    collection_date_hist(df_met, svfld / "collection_date_hist.png", dpi=300)

    # metadata tree
    metadata_tree(
        df_met,
        df_all,
        tree,
        svfld / "metadata_tree.png",
        svfld / "metadata_tree_legend.png",
        dpi=300,
    )
