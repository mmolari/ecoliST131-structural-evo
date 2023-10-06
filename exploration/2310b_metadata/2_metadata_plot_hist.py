# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib as mpl
from Bio import Phylo

# %%
tree_file = "../../results/ST131_full/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
tree.ladderize()

# dff = "../../results/ST131/metadata.csv"
dff = "../../results/ST131_full/metadata.csv"
df1 = pd.read_csv(dff, index_col=0)
dff = "../../results/ST131_full/assembly_qc/alleles_summary.csv"
df2 = pd.read_csv(dff, index_col=0)
df = pd.concat([df1, df2], axis=1)
# %%

# histogram
srs = df["geo_loc_name"]
fig, ax = plt.subplots(figsize=(8, 6))
srs.value_counts().plot(kind="bar", ax=ax)
# write number without value
n_nan = srs.isna().sum()
ax.text(0.8, 0.9, f"# NaN = {n_nan}", transform=ax.transAxes)
sns.despine()
plt.tight_layout()
plt.show()


srs = df["continent"]
fig, ax = plt.subplots(figsize=(8, 6))
srs.value_counts().plot(kind="bar", ax=ax)
# write number without value
n_nan = srs.isna().sum()
ax.text(0.8, 0.9, f"# NaN = {n_nan}", transform=ax.transAxes)
sns.despine()
plt.tight_layout()
plt.show()

# %%

srs = df["year"]
fig, ax = plt.subplots(figsize=(8, 6))
srs.value_counts().plot(kind="bar", ax=ax)
# write number without value
n_nan = srs.isna().sum()
ax.text(0.8, 0.9, f"# NaN = {n_nan}", transform=ax.transAxes)
sns.despine()
plt.tight_layout()
plt.show()
# %%

srs = df["host"]
fig, ax = plt.subplots(figsize=(8, 6))
srs.value_counts().plot(kind="bar", ax=ax)
# write number without value
n_nan = srs.isna().sum()
ax.text(0.8, 0.9, f"# NaN = {n_nan}", transform=ax.transAxes)
sns.despine()
plt.tight_layout()
plt.show()
# %%

srs = df["isolation_source"]
fig, ax = plt.subplots(figsize=(8, 6))
srs.value_counts().plot(kind="bar", ax=ax)
# write number without value
n_nan = srs.isna().sum()
ax.text(0.8, 0.9, f"# NaN = {n_nan}", transform=ax.transAxes)
sns.despine()
plt.tight_layout()
plt.show()
# %%

# country histogram
ctrs = df[["geo_loc_name", "continent"]].value_counts()

# assign color to each continent
color_dict = {
    "Africa": "tab:purple",
    "Asia": "tab:red",
    "Europe": "tab:blue",
    "North America": "tab:green",
    "Oceania": "tab:orange",
    "South America": "tab:brown",
}

continent_ct = df["continent"].value_counts()
# create legend
legend_elements = []
for k, v in color_dict.items():
    legend_elements.append(
        mpl.lines.Line2D(
            [0], [0], color=v, label=k + f" (n={continent_ct[k]})", marker="s", ls=""
        )
    )

fig, ax = plt.subplots(1, 1, figsize=(10, 5))
# draw colored barplot
xticks = []
xticklabels = []
for idx, row in ctrs.iteritems():
    ctr, cont = idx
    ax.bar(ctr, row, color=color_dict[cont])
# set labels
ax.set_xlabel("Country")
ax.set_ylabel("Number of isolates")
# rotate xticks
plt.setp(ax.get_xticklabels(), rotation=90, ha="center")
ax.legend(handles=legend_elements, loc="upper right", title="Continent")
# %%

# plot collection date


def collection_date_hist(df):
    yrs = df["year"]
    M, m = yrs.max(), yrs.min()
    bins = np.arange(m, M + 3) - 0.5
    fig, ax = plt.subplots(
        1,
        1,
    )
    ax.hist(yrs, bins=bins, edgecolor="black", color="lightgrey")
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of isolates")
    sns.despine()
    plt.tight_layout()
    plt.show()


collection_date_hist(df)

# %%
continent_color = {
    "Africa": "tab:purple",
    "Asia": "tab:red",
    "Europe": "tab:blue",
    "North America": "tab:green",
    "Oceania": "tab:orange",
    "South America": "tab:brown",
}


def __mark(ax, x, y, c, **kwargs):
    ax.scatter(x, y, marker="s", edgecolor="lightgray", color=c, **kwargs)


def create_legend(color_dictionaries):
    L = len(color_dictionaries)
    fig, axs = plt.subplots(1, L, figsize=(L * 2, 2))
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
    plt.show()


def metadata_tree(df, alleles, tree):
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
    plt.show()

    color_dicts = {
        "country": country_color,
        "source": isolation_color,
        "year": year_color,
    } | allele_colors
    create_legend(color_dicts)


alleles = df2.columns.to_list()
metadata_tree(df, alleles, tree)

# %%
