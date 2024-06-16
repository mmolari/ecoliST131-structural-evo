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


def parse_input_data():

    # tree
    tree_file = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
    tree = Phylo.read(tree_file, "newick")
    strains = [x.name for x in tree.get_terminals()]
    strain_y = {s: i + 1 for i, s in enumerate(strains)}

    # metadata
    metadata_file = "../../results/ST131_ABC/metadata.csv"
    mdf = pd.read_csv(metadata_file, index_col=0)

    # MLST chromosome
    alleles_file = "../../results/ST131_ABC/assembly_qc/alleles_summary.csv"
    adf = pd.read_csv(alleles_file, index_col=0)

    # plasmid MLST
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

    # resistance
    resistance_file = (
        # "../../results/ST131_ABC/plasmids/resistance/ncbi_summary_chromosome.csv"
        "../../results/ST131_ABC/plasmids/resistance/card_summary_chromosome.csv"
    )
    rdf_pld = pd.read_csv(resistance_file, index_col=0)
    # resistance_file = "../../results/ST131_ABC/resistance/ncbi_summary.txt"
    resistance_file = "../../results/ST131_ABC/resistance/card_summary.txt"
    rdf_chr = pd.read_csv(resistance_file, index_col=0, sep="\t")
    rdf_chr.index = rdf_chr.index.str.split("/").str[-1].str.removesuffix(".tab")
    rdf_chr = rdf_chr.drop(columns=["NUM_FOUND"])

    def convert(x):
        if x == ".":
            return 0
        elif isinstance(x, str):
            x = np.array([float(i) for i in x.split(";")])
            x = np.sum(x > 80)
            return x
        else:
            return int(x > 80)

    rdf_chr = rdf_chr.map(convert)

    return tree, strains, strain_y, mdf, adf, pldf, plasmid_data, rdf_chr, rdf_pld


tree, strains, strain_y, mdf, adf, pldf, plasmid_data, rdf_chr, rdf_pld = (
    parse_input_data()
)


# %%
def draw_legend(colors, svname, title=None):
    figl, axl = plt.subplots(figsize=(3, 3))
    for cont, col in colors.items():
        axl.plot([], [], "s", color=col, label=cont)
    axl.legend(title=title, loc="center")
    axl.set_axis_off()
    figl.tight_layout()
    figl.savefig(str(svname) + ".pdf", bbox_inches="tight")
    figl.savefig(str(svname) + ".svg", bbox_inches="tight")
    plt.close(figl)


def draw_continent(ax, x, mdf, strains):
    colors = {
        cont: col
        for cont, col in zip(
            mdf["continent"].dropna().unique(), sns.color_palette("Set2")
        )
    }
    for s in strains:
        cont = mdf.loc[s, "continent"]
        if cont is np.nan:
            continue
        y = strain_y[s]
        col = colors[cont]
        ax.plot(
            x,
            y,
            "s",
            color=col,
            # markeredgecolor="black",
            markersize=4,
            markeredgewidth=0.5,
        )
    # legend
    draw_legend(colors, fig_fld / "leg_continent", title="Continent")


def draw_year(ax, x, mdf, strains):
    year_cat = pd.cut(
        mdf["year"],
        [2000, 2010, 2020, 2030],
        labels=["2000-2009", "2010-2019", "2020-2022"],
        right=False,
    )
    colors = {
        "2000-2009": sns.color_palette("Blues")[1],
        "2010-2019": sns.color_palette("Blues")[3],
        "2020-2022": sns.color_palette("Blues")[5],
    }
    for s in strains:
        year = year_cat.loc[s]
        if year is np.nan:
            continue
        y = strain_y[s]
        col = colors[year]
        ax.plot(
            x,
            y,
            "s",
            color=col,
            markeredgecolor="none",
            markersize=4,
            markeredgewidth=0.5,
        )
    # legend
    draw_legend(colors, fig_fld / "leg_year", title="Year")


def draw_allele(ax, x, adf, strains):
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

    colors = {str(v).split(".")[0]: c for v, c in colors.items()}

    draw_legend(colors, fig_fld / f"leg_all_{al}", title=al.removesuffix("_eb"))


def draw_chr_resistance(ax, x, strains, rdf, gene, color):
    rs = rdf[gene]
    for s in strains:
        r = rs.loc[s]
        if (r is np.nan) or (r < 1):
            continue
        y = strain_y[s]
        ax.plot(
            x,
            y,
            "o",
            markerfacecolor="none",
            markeredgecolor=color,
            markersize=4,
            markeredgewidth=0.5,
        )


def draw_pld_resistance(ax, x, strains, rdf, gene, color):
    rs = rdf[gene]
    for s in strains:
        if not s in plasmid_data:
            continue
        r = rs.loc[s]
        if (r is np.nan) or (r < 1):
            continue
        y = strain_y[s]
        ax.plot(
            x,
            y,
            "x",
            markerfacecolor="none",
            markeredgecolor=color,
            markersize=4,
            markeredgewidth=0.5,
        )


def draw_resistance_leg(svname, c_chrom, c_plasmid):
    figl, axl = plt.subplots(figsize=(3, 3))
    axl.plot(
        [],
        [],
        "o",
        markerfacecolor="none",
        markeredgecolor=c_chrom,
        markersize=4,
        markeredgewidth=0.5,
        label="on chromosome",
    )
    axl.plot(
        [],
        [],
        "x",
        markerfacecolor="none",
        markeredgecolor=c_plasmid,
        markersize=4,
        markeredgewidth=0.5,
        label="on plasmid",
    )
    axl.legend(title="Resistance", loc="center")
    axl.set_axis_off()
    figl.tight_layout()
    figl.savefig(str(svname) + ".pdf", bbox_inches="tight")
    figl.savefig(str(svname) + ".svg", bbox_inches="tight")
    plt.close(figl)


# %%
fig, axs = plt.subplots(
    1, 2, figsize=(6, 10), sharey=True, gridspec_kw={"width_ratios": [1, 1.2]}
)

# draw tree
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
xticks, xlabels = [], []
# continent
x = 0
draw_continent(ax, x, mdf, strains)
xticks.append(x)
xlabels.append("Continent")

# years
x = 1
draw_year(ax, x, mdf, strains)
xticks.append(x)
xlabels.append("Year")

# alleles
x = 3
for al in ["fimH_eb", "gyrA_eb", "parC_eb"]:
    draw_allele(ax, x, adf, strains)
    xticks.append(x)
    xlabels.append(al.removesuffix("_eb"))
    x += 1


# plasmid MLST
x = 7
values = pldf["typ"].dropna().value_counts().index[:6].to_list()
# values = [
#     "F29:A-:B10",
#     "F29:A-:B-",
#     "F2:A1:B-",
#     "F2:A-:B-",
#     "F1:A2:B20",
#     "F1:A2:B-",
#     "F36:A4:B-",
# ]

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
x = 14
c_chr = "#2a9d8f"
c_plm = "#f4a261"
values = ["blaCTX-M-14", "blaCTX-M-15", "blaCTX-M-27", "blaOXA-1", "aac(6')-Ib-D181Y"]
values = ["CTX-M-14", "CTX-M-15", "CTX-M-27", "OXA-1", "AAC(6')-Ib-cr"]
for v in values:
    if v in rdf_chr.columns:
        draw_chr_resistance(ax, x, strains, rdf_chr, v, c_chr)
    if v in rdf_pld.columns:
        draw_pld_resistance(ax, x, strains, rdf_pld, v, c_plm)
    xticks.append(x)
    xlabels.append(v)
    x += 1
draw_resistance_leg(svname=fig_fld / "leg_res", c_chrom=c_chr, c_plasmid=c_plm)


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
fig.savefig(fig_fld / "tree.pdf", bbox_inches="tight")
fig.savefig(fig_fld / "tree.svg", bbox_inches="tight")
plt.show()

# %%

# %%

tree_len = tree.total_branch_length()
print(f"n. mutations per variant = ", tree_len * 2427416 / 1936)
# %%
