# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import pathlib
import json
from Bio import Phylo

fig_fld = pathlib.Path("figs/f04")
fig_fld.mkdir(exist_ok=True, parents=True)

data_fld = pathlib.Path("data/f04")
data_fld.mkdir(exist_ok=True, parents=True)

dset = "ST131_ABC"
fld = pathlib.Path(f"../../results/{dset}")

# load annotations
fld_data = fld / "annotations/junct_pos_asm20-100-5"
is_df = pd.read_csv(fld_data / "ISEScan_real.csv", index_col=0)

# load info
is_info = pd.read_csv(fld / "annotations/isescan/is_summary.tsv", sep="\t")
is_info["id"] = (
    is_info["seqID"] + "|" + is_info["family"] + "|" + is_info.index.astype(str)
)
is_info = is_info.set_index("id", verify_integrity=True)


# load joint coordinates dictionary
pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
with open(pos_file) as f:
    jp = json.load(f)

# load junction info dataframes
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_j = pd.read_csv(df_file, index_col=0)
df_j["delta_len"] = df_j["max_length"] - df_j["min_length"]

# load tree
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]

# genome lengths
iso_L_file = fld / "pangraph/genome_lengths.csv"
iso_L = pd.read_csv(iso_L_file, index_col=0)["length"].to_dict()

# %%

# add number of IS
df_j["n_IS"] = is_df["junction"].value_counts()
df_j["n_IS"] = df_j["n_IS"].fillna(0)
df_j["n_IS"] = df_j["n_IS"].astype(int)

mask = df_j["n_IS"] > 0
mask &= df_j["n_categories"] == 2
mask &= df_j["n_iso"] == 222
sdf = df_j[mask]
Js = sdf.index.to_list()
sdf

# %%
summary_df = pd.read_csv(data_fld / "summary.csv", index_col=0)
summary_df["n_products"] = summary_df["n_products"].fillna(0)
summary_df

# %%
cols = [
    "junction",
    "n_empty",
    "n_full",
    "n_IS",
    "n_IS_isolates",
    "nonempty_is_IS",
    "n_hits_gene",
    "n_hits_CDS",
    "n_products",
    "annotation_consistency",
    "annotations",
]

summary_df[["nonempty_is_IS", "n_products", "annotation_consistency"]].value_counts(
    dropna=False
)

# %%
mask = summary_df["nonempty_is_IS"] & (summary_df["n_products"] > 1)
summary_df[mask]

# %%

mask = summary_df["nonempty_is_IS"] & summary_df["n_products"].isna()
summary_df[mask]["n_hits_gene"].value_counts()

# %%

# %%

# flow plot

tot = {
    "total": len(summary_df),
    "color": "slateblue",
    "title": "IS-coldspots",
}

consistent = summary_df["nonempty_is_IS"].value_counts().to_dict()
C, NC = consistent[True], consistent[False]
layer_1 = {
    "values": [C, NC],
    "labels": ["consistent", "inconsistent"],
    "colors": ["C0", "gray"],
    "title": "IS/non-empty consistency",
    "from": [[(0, C)], [(0, NC)]],
}
products = (
    summary_df[["nonempty_is_IS", "n_products"]].value_counts(dropna=False).to_dict()
)
C0 = products[(True, 0)]
C1 = products[(True, 1)]
C2 = products[(True, 2)]
N0 = products[(False, 0)]

layer_2 = {
    "values": [C1 + C2, C0],
    "labels": ["yes", "no"],
    "colors": ["C2", "C3"],
    "title": "has annotations",
    "from": [[(0, C1 + C2)], [(0, C0)]],
}

layer_3 = {
    "values": [C1, C2],
    "labels": ["n=1", "n=2"],
    "colors": ["limegreen", "darkgreen"],
    "title": "n. products",
    "from": [[(0, C1)], [(0, C2)]],
}

layers = [layer_1, layer_2, layer_3]


def draw_layer(layer, x, w, tot, connect, ax):

    V = layer["values"]
    Y = []
    y = 0
    y_low_left = {}
    for i, v in enumerate(V):
        c = layer["colors"][i]
        lab = layer["labels"][i]
        ax.bar(x, v, w, y, color=c, linewidth=1, edgecolor="k")
        ax.text(
            x - 0.5 * (1 - w),
            y + v / 2,
            # f"{v} ({v/tot:.1%})\n{lab}",
            f"n={v}\n{lab}",
            ha="center",
            va="center",
        )
        Y.append((y, y + v))
        F = layer["from"][i]
        y_low_right = y
        for idx, h in F:
            X = [x - 1 + w * 0.5, x - 0.5 * w]
            if not idx in y_low_left:
                y_low_left[idx] = connect[idx][0]
            y_high_right = y_low_right + h
            y_high_left = y_low_left[idx] + h
            ax.fill_between(
                X,
                [y_low_left[idx], y_low_right],
                [y_high_left, y_high_right],
                color=c,
                alpha=0.2,
                linewidth=0,
            )
            y_low_right += h
            y_low_left[idx] += h

        y += v + 10
    return Y


def draw_flowplot(tot, layers, w=0.2, figsize=(8, 5)):
    fig, ax = plt.subplots(figsize=figsize)

    x = 0
    T = tot["total"]
    c = tot["color"]
    ax.bar(x, T, w, 0, color=c, label="total", linewidth=1, edgecolor="k")
    ax.text(x - 0.5 * (1 - w), T / 2, f"n={T}", ha="center", va="center")
    x += 1
    connect = [(0, T)]

    for layer in layers:
        connect = draw_layer(layer, x, w, T, connect, ax)
        x += 1

    xticks = [tot["title"]] + [layer["title"] for layer in layers]
    ax.set_xticks(range(x))
    ax.set_xticklabels(xticks)
    ax.set_xlim(left=-1 + w)
    return fig, ax


fig, ax = draw_flowplot(tot, layers)
ax.set_ylabel("n. junctions")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "flowplot.png")
plt.show()
# %%
mask = summary_df["nonempty_is_IS"] & (summary_df["n_products"] == 1)
summary_df[mask]["annotations"].to_csv(data_fld / "broken_products.csv")
# %%


# flow plot

# backbone coldspots
mask = df_j["n_categories"] == 2
mask &= df_j["n_iso"] == 222


tot = {
    "total": mask.sum(),
    "color": "slateblue",
    "title": "binary junctions",
}

Iass = len(summary_df)
NIass = mask.sum() - Iass
layer_0 = {
    "values": [Iass, NIass],
    "labels": ["IS", "no IS"],
    "colors": ["C0", "gray"],
    "title": "IS-association",
    "from": [[(0, Iass)], [(0, NIass)]],
}

consistent = summary_df["nonempty_is_IS"].value_counts().to_dict()
C, NC = consistent[True], consistent[False]
layer_1 = {
    "values": [C, NC],
    "labels": ["consistent", "inconsistent"],
    "colors": ["C0", "gray"],
    "title": "IS/non-empty consistency",
    "from": [[(0, C)], [(0, NC)]],
}
products = (
    summary_df[["nonempty_is_IS", "n_products"]].value_counts(dropna=False).to_dict()
)
C0 = products[(True, 0)]
C1 = products[(True, 1)]
C2 = products[(True, 2)]
N0 = products[(False, 0)]

layer_2 = {
    "values": [C1 + C2, C0],
    "labels": ["yes", "no"],
    "colors": ["C2", "C3"],
    "title": "has annotations",
    "from": [[(0, C1 + C2)], [(0, C0)]],
}

layer_3 = {
    "values": [C1, C2],
    "labels": ["1 product", "2 products"],
    "colors": ["limegreen", "darkgreen"],
    "title": "n. products",
    "from": [[(0, C1)], [(0, C2)]],
}

layers = [layer_0, layer_1, layer_2, layer_3]

# %%


fig, ax = draw_flowplot(tot, layers, figsize=(9, 5))
ax.set_ylabel("n. junctions")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "flowplot_coldspots.png")
plt.savefig(fig_fld / "flowplot_coldspots.pdf")
plt.show()

# %%
