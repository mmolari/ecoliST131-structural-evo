# %%
from Bio import Phylo
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pathlib
import json

fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(parents=True, exist_ok=True)

# %%
tree_file = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]
strain_y = {s: i + 1 for i, s in enumerate(strains)}

metadata_file = "../../results/ST131_ABC/metadata.csv"
mdf = pd.read_csv(metadata_file, index_col=0)

alleles_file = "../../results/ST131_ABC/assembly_qc/alleles_summary.csv"
adf = pd.read_csv(alleles_file, index_col=0)

plasmid_mlst_file = "../../results/ST131_ABC/plasmids/alleles_summary.csv"
pldf = pd.read_csv(plasmid_mlst_file, index_col=0)


def parse_mlst_df(mlst_df_file):
    df = pd.read_csv(mlst_df_file, sep=",", index_col=0, dtype=str)
    df.fillna("-", inplace=True)
    df["typ"] = "F" + df["FII"] + ":A" + df["FIA"] + ":B" + df["FIB"]

    # remove untyped
    mask = df["typ"] == "F-:A-:B-"
    df = df[~mask]
    return df


pldf = parse_mlst_df(plasmid_mlst_file)

plasmid_json_file = "../../config/datasets/ST131_ABC/plasmids.json"
with open(plasmid_json_file, "r") as f:
    plasmid_data = json.load(f)

resistance_file = (
    "../../results/ST131_ABC/plasmids/resistance/ncbi_summary_chromosome.csv"
)
rdf_chr = pd.read_csv(resistance_file, index_col=0)
resistance_file = "../../results/ST131_ABC/plasmids/resistance/ncbi_summary_plasmid.txt"
rdf_pld = pd.read_csv(resistance_file, index_col=0, sep="\t")

# %%
fig, axs = plt.subplots(
    1, 2, figsize=(8, 10), sharey=True, gridspec_kw={"width_ratios": [2, 1]}
)

ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
right = ax.get_xlim()[1]
ax.set_xlabel("")
ax.set_ylabel("")

for i, strain in enumerate(strains):
    d = tree.distance(strain)
    ax.plot([d, right], [i + 1, i + 1], ":", color="gray", lw=0.5)

ax.set_xlim(right=right * 0.9)

ax = axs[1]
xticks = []
xlabels = []
# continent
x = 0
colors = {
    cont: col
    for cont, col in zip(mdf["continent"].dropna().unique(), sns.color_palette("Set2"))
}
for s in strains:
    cont = mdf.loc[s, "continent"]
    if cont is np.nan:
        continue
    y = strain_y[s]
    col = colors[cont]
    ax.plot(
        x, y, "s", color=col, markeredgecolor="black", markersize=4, markeredgewidth=0.5
    )
xticks.append(x)
xlabels.append("Continent")
# legend
figl, axl = plt.subplots()
for cont, col in colors.items():
    axl.plot([], [], "s", color=col, label=cont)
axl.legend(title="Continent", loc="center")
axl.set_axis_off()


# years
x = 1
year_cat = pd.cut(
    mdf["year"],
    [2000, 2010, 2020, 2030],
    labels=["2000-2010", "2010-2020", ">2020"],
)
colors = {
    "2000-2010": sns.color_palette("Blues")[0],
    "2010-2020": sns.color_palette("Blues")[2],
    ">2020": sns.color_palette("Blues")[5],
}
for s in strains:
    year = year_cat.loc[s]
    if year is np.nan:
        continue
    y = strain_y[s]
    col = colors[year]
    ax.plot(
        x, y, "s", color=col, markeredgecolor="black", markersize=4, markeredgewidth=0.5
    )
xticks.append(x)
xlabels.append("Year")
# legend
figl, axl = plt.subplots()
for year, col in colors.items():
    axl.plot([], [], "s", color=col, label=year)
axl.legend(title="Year", loc="center")
axl.set_axis_off()

ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, rotation=90)

# alleles
x = 3
for al in ["fimH_eb", "gyrA_eb", "parC_eb"]:

    values = adf[al].dropna().value_counts().index[:3]
    colors = sns.color_palette("tab10", n_colors=len(values))
    colors = {v: c for v, c in zip(values, colors)}
    colors["other"] = "lightgray"

    for s in strains:
        a = adf.loc[s, al]
        if a is np.nan:
            continue
        c = colors[a] if a in colors else colors["other"]
        y = strain_y[s]
        ax.plot(x, y, "s", color=c, markersize=4, markeredgewidth=0.5)
    xticks.append(x)
    xlabels.append(al)
    figl, axl = plt.subplots()
    for value, col in colors.items():
        try:
            v = int(float(value))
        except:
            v = value
        axl.plot([], [], "s", color=col, label=v)
    axl.legend(title=al, loc="center")

    x += 1


# plasmid MLST
x = 7
values = pldf["typ"].dropna().value_counts().index[:6]

for v in values:
    for s in strains:
        if not s in plasmid_data:
            continue
        plasmids = plasmid_data[s]
        is_there = False
        for p in plasmids:
            if p in pldf.index:
                if pldf.loc[p]["typ"] == v:
                    is_there = True
                    break
        if is_there:
            ax.plot(
                x, strain_y[s], "s", color="black", markersize=4, markeredgewidth=0.5
            )
    xticks.append(x)
    xlabels.append(v)
    x += 1

# resistance
x = 13
values = ["blaCTX-M-14", "blaCTX-M-15", "blaCTX-M-27"]


# final setup
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, rotation=90)

xl, xr = ax.get_xlim()
for s in strains:
    y = strain_y[s]
    ax.plot([xl, xr], [y, y], ":", color="gray", lw=0.5, zorder=-1)
ax.set_xlim(left=xl, right=xr)
fig.subplots_adjust(wspace=0.0)
sns.despine(fig)
# set yticks invisible
ax.spines["left"].set_visible(False)
ax.yaxis.set_visible(False)


plt.tight_layout()
plt.show()

# %%

# %%
