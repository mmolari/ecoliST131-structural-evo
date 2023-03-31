# %%
import pathlib
import numpy as np
import pandas as pd
from Bio import Phylo
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

fig_fld = pathlib.Path("figs")
fig_fld.mkdir(exist_ok=True)


def svfig(svname):
    for suffix in ["png", "pdf", "svg"]:
        plt.savefig(fig_fld / f"{svname}.{suffix}", dpi=300)


# %%

metadata_fix = {
    "collection_date": {"2015-01/2016-07": "2015"},
    "host": {
        "pantropical spotted dolphin": "dolphin",
        "Homo sapiens; 55 year old": "human",
        "Pig": "pig",
        "Homo sapiens": "human",
    },
    "geo_loc_name": {
        "USA: Minnesota": "USA",
        "USA: New York": "USA",
        "USA: Minneapolis": "USA",
        "USA: Burlington": "USA",
        "unknown": pd.NA,
        "Canada: Winnipeg": "Canada",
        "Brazil: Botucatu, Sao Paulo State": "Brazil",
        "USA:Utah": "USA",
        "USA:Houston": "USA",
        "United Kingdom:Manchester": "United Kingdom",
        "Australia:Brisbane": "Australia",
        "Czech Republic: Brno": "Czech Republic",
        "USA: San Diego": "USA",
        "USA: Houston": "USA",
        "China: Hainan Province": "China",
        "USA: Pittsburgh, Pennsylvania": "USA",
    },
}

# %%
df = pd.read_csv("../../metadata/ST131/metadata-raw.csv", index_col=0)
for k, d in metadata_fix.items():
    df[k] = df[k].replace(d)

mask = df["isolation_source"].isin(
    ["wastewater", "wastewater treatment plant effluent"]
)
df.loc[mask, "host"] = "wastewater"

df["collection_date"] = pd.to_datetime(df["collection_date"])
# df["year"] = df["collection_date"].map(
#     lambda x: str(x.year) if not pd.isna(x) else np.nan
# )
# order_yr = np.sort(df["year"].value_counts().index)

df["year"] = df["collection_date"].map(
    lambda x: int(x.year) if not pd.isna(x) else np.nan
)
df["year"] = pd.cut(
    df["year"],
    bins=[2002, 2008, 2012, 2016, 2020],
    include_lowest=True,
    labels=["2002-2008", "2009-2012", "2013-2016", "2017-2020"],
)

order_geo = df["geo_loc_name"].value_counts().index
df["geo_loc_name"] = pd.Categorical(df["geo_loc_name"], order_geo, ordered=True)

order_host = df["host"].value_counts().index
df["host"] = pd.Categorical(df["host"], order_host, ordered=True)

# %%


def addna(ax, series):
    n = series.isna().sum()
    ax.text(0.7, 0.1, f"N.A. = {n}", transform=ax.transAxes)


fig, axs = plt.subplots(1, 3, figsize=(11, 4))

ax = axs[0]
sns.histplot(data=df, y="host", ax=ax)
addna(ax, df["host"])
ax.set_ylabel("host")

ax = axs[1]
sns.histplot(data=df, y="collection_date", ax=ax)
addna(ax, df["collection_date"])
ax.set_ylabel("collection date")


ax = axs[2]
sns.histplot(data=df, y="geo_loc_name", ax=ax)
addna(ax, df["geo_loc_name"])
ax.set_ylabel("location")

sns.despine()
plt.tight_layout()
svfig("metadata_hist")
plt.show()
# %%
tree = Phylo.read(
    "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk", format="newick"
)
tree.ladderize()
str_order = [x.name for x in tree.get_terminals()]
N = len(str_order)

# %%
# assign colors
c_dict = {}

cmap = plt.get_cmap("Pastel1")
cat = df["host"].dtype.categories
c_dict["host"] = {v: cmap(n) for n, v in enumerate(cat)}
c_dict["host"][np.nan] = "white"

cmap = plt.get_cmap("viridis")
cat = df["year"].dtype.categories
c_dict["year"] = {v: cmap(n / (len(cat) - 1)) for n, v in enumerate(cat)}
c_dict["year"][np.nan] = "white"

cmap = plt.get_cmap("tab20")
cat = df["geo_loc_name"].dtype.categories
c_dict["geo_loc_name"] = {v: cmap((n*2 % 20) + (n > 9)) for n, v in enumerate(cat)}
c_dict["geo_loc_name"][np.nan] = "white"
# %%


fig, axs = plt.subplots(
    1,
    2,
    figsize=(8, 8),
    # sharey=True,
    gridspec_kw={"width_ratios": [3, 1]},
)

ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: None)
ax.grid(True, alpha=0.3)
xl, xr = ax.get_xlim()
ddx = (xr - xl) * 0.05
for i, k in enumerate(["host", "year", "geo_loc_name"]):
    D = df[k].to_dict()
    x = xr + ddx * np.ones(N) * (i - 2)
    y = np.arange(1, N + 1)
    c = [c_dict[k][D[s]] for s in str_order]
    ax.scatter(x, y, color=c, marker="o")

ax.set_xlim(right=xr + ddx)

for ax in axs:
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    for k in ax.spines:
        ax.spines[k].set_visible(False)

ax = axs[1]
ax.set_yticks([])

# legends
E = []
for i, k in enumerate(["host", "year", "geo_loc_name"]):
    elm = [
        mpl.lines.Line2D([], [], marker="o", ls="", color=c, label=l)
        for l, c in c_dict[k].items()
        if isinstance(l, str)
    ]
    E.append(elm)
l0 = ax.legend(handles=E[0], title="host", loc="upper left")
l1 = ax.legend(handles=E[1], title="year", loc="center left")
l2 = ax.legend(handles=E[2], title="location", loc="lower left")
ax.add_artist(l0)
ax.add_artist(l1)
plt.tight_layout()
svfig("metadata_tree")
plt.show()

# %%
