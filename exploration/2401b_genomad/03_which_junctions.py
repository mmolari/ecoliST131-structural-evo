# %%
import pathlib
import json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import location_utils as lu
from Bio import SeqIO

fig_fld = pathlib.Path("figs/f03")
fig_fld.mkdir(exist_ok=True, parents=True)

data_fld = pathlib.Path("data")
data_fld.mkdir(exist_ok=True, parents=True)

dset = "ST131_ABC"


def preprocess(df):
    df["iso"] = df["seq_name"].str.split("|").str[0]
    df["start"] = df["coordinates"].str.split("-").str[0].astype(int)
    df["end"] = df["coordinates"].str.split("-").str[1].astype(int)
    df = df.drop(columns=["fdr", "taxonomy", "topology", "genetic_code", "coordinates"])
    return df


# load prophage dataframe
fld = pathlib.Path(f"../../results/{dset}")
fname = fld / "annotations/genomad/prophage_summary.tsv"
df_ph = pd.read_csv(fname, sep="\t")
df_ph = preprocess(df_ph)
df_ph.set_index("seq_name", inplace=True)

# load junction info dataframes
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_j = pd.read_csv(df_file, index_col=0)
df_j["delta_len"] = df_j["max_length"] - df_j["min_length"]

# load edge length dataframe
df_edge_len = fld / "backbone_joints/asm20-100-5/edge_len.csv"
df_el = pd.read_csv(df_edge_len, index_col=0)

# load joint coordinates dictionary
pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
with open(pos_file) as f:
    jp = json.load(f)

# isolates genome length
strains = df_el.index.to_list()
iso_L = {}
for iso in strains:
    fa_fname = fld / f"../../data/fa/{iso}.fa"
    with open(fa_fname) as f:
        iso_L[iso] = len(SeqIO.read(f, "fasta"))

# %%
pp_to_j = []
rand_pp_to_j = []

for pp, row in df_ph.iterrows():
    pb, pe = row["start"], row["end"]
    iso = row["iso"]
    # locate on joints
    for j, pos_d in jp.items():
        if not iso in pos_d:
            continue
        cb, ab, ae, ce, strand = pos_d[iso]
        in_core = lu.within(cb, ce, pb, pe)
        if in_core == "no":
            continue
        in_acc = lu.within(ab, ae, pb, pe)
        if in_acc == "no":
            continue
        pp_to_j.append(
            {
                "phage": pp,
                "junction": j,
                "iso": iso,
                "in_core": in_core,
                "in_acc": in_acc,
                "pb": pb,
                "pe": pe,
                "jcb": cb,
                "jce": ce,
                "jab": ab,
                "jae": ae,
                "strand": strand,
            }
        )
    # locate on random joints
    gL = iso_L[iso]
    re_place = True
    while re_place:
        delta_L = np.random.randint(gL)
        pbn, pen = (pb + delta_L) % gL, (pe + delta_L) % gL
        if pbn < pen:
            re_place = False
    for j, pos_d in jp.items():
        if not iso in pos_d:
            continue
        cb, ab, ae, ce, strand = pos_d[iso]
        in_core = lu.within(cb, ce, pbn, pen)
        if in_core == "no":
            continue
        in_acc = lu.within(ab, ae, pbn, pen)
        if in_acc == "no":
            continue
        rand_pp_to_j.append(
            {
                "phage": pp,
                "junction": j,
                "iso": iso,
                "in_core": in_core,
                "in_acc": in_acc,
                "pb": pb,
                "pe": pe,
                "jcb": cb,
                "jce": ce,
                "jab": ab,
                "jae": ae,
                "strand": strand,
            }
        )


pp_to_j = pd.DataFrame(pp_to_j)
pp_to_j.to_csv(data_fld / "phages_to_joints.csv")
rand_pp_to_j = pd.DataFrame(rand_pp_to_j)
pp_to_j
# %%

df_ph["assigned"] = pp_to_j["phage"].value_counts()
df_ph["assigned"].fillna(0, inplace=True)
df_ph["assigned_rand"] = rand_pp_to_j["phage"].value_counts()
df_ph["assigned_rand"].fillna(0, inplace=True)

sns.histplot(
    data=df_ph,
    x="assigned",
    discrete=True,
    element="step",
    label="geNomad annotations",
)
sns.histplot(
    data=df_ph,
    x="assigned_rand",
    discrete=True,
    element="step",
    label="shuffled annotations",
    alpha=0.5,
    color="gray",
)
plt.legend()
plt.xlabel("number of joints")
plt.ylabel("number of prophage annotations")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "joints_per_phage.png", dpi=300)
plt.show()

# %%
df_j["n_prophages"] = pp_to_j["junction"].value_counts()
df_j["n_prophages"].fillna(0, inplace=True)
df_j["has_prophage"] = df_j["n_prophages"] > 0
df_j["prophage_freq"] = (
    pp_to_j.groupby("junction").apply(lambda x: len(np.unique(x["iso"])))
    / df_j["n_iso"]
)
df_j["prophage_freq"].fillna(0, inplace=True)

df_j["n_prophages_rand"] = rand_pp_to_j["junction"].value_counts()
df_j["n_prophages_rand"].fillna(0, inplace=True)
df_j["has_prophage_rand"] = df_j["n_prophages_rand"] > 0
df_j["prophage_freq_rand"] = (
    rand_pp_to_j.groupby("junction").apply(lambda x: len(np.unique(x["iso"])))
    / df_j["n_iso"]
)
df_j["prophage_freq_rand"].fillna(0, inplace=True)


sns.histplot(
    data=df_j, x="n_prophages", discrete=True, element="step", label="geNomad joints"
)
sns.histplot(
    data=df_j,
    x="n_prophages_rand",
    discrete=True,
    element="step",
    label="shuffled joints",
    alpha=0.5,
    color="gray",
)
plt.xscale("symlog")
plt.legend()
plt.xlabel("number of prophage annotations")
plt.ylabel("number of joints")
plt.yscale("log")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "phage_per_joint.png", dpi=300)
plt.show()


# %%
sns.scatterplot(
    data=df_j,
    x="nonempty_freq",
    y="prophage_freq",
    hue="has_prophage",
    alpha=0.3,
)
plt.plot([0, 1], [0, 1], color="gray", linestyle="--", zorder=-3, alpha=0.2)
plt.xlabel("nonempty frequency")
plt.ylabel("prophage frequency")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "phage_vs_nonempty_freq.png", dpi=300)
plt.show()

# sns.scatterplot(
#     data=df_j,
#     x="nonempty_freq",
#     y="prophage_freq",
#     alpha=0.3,
#     label="geNomad annotations",
# )
# sns.scatterplot(
#     data=df_j,
#     x="nonempty_freq",
#     y="prophage_freq_rand",
#     alpha=0.3,
#     label="shuffled annotations",
# )
# plt.plot([0, 1], [0, 1], color="gray", linestyle="--", zorder=-3, alpha=0.2)
# plt.xlabel("nonempty frequency")
# plt.ylabel("prophage frequency")
# plt.legend()
# sns.despine()
# plt.tight_layout()
# plt.savefig(fig_fld / "phage_vs_nonempty_freq_rand.png", dpi=300)
# plt.show()

# %%

# joint plot
fig, axs = plt.subplots(
    2,
    2,
    figsize=(6, 6),
    sharex="col",
    sharey="row",
    gridspec_kw={"width_ratios": [1, 0.2], "height_ratios": [0.2, 1]},
)

ax = axs[1, 0]
sns.scatterplot(
    data=df_j,
    x="delta_len",
    y="n_categories",
    hue="has_prophage",
    alpha=0.2,
    ax=ax,
)
ax.set_xlabel("delta length")
ax.set_ylabel("n categories")

sns.histplot(
    data=df_j,
    x="delta_len",
    hue="has_prophage",
    alpha=0.2,
    ax=axs[0, 0],
    legend=False,
    log_scale=True,
    bins=25,
    element="step",
    common_bins=True,
    # common_norm=False,
    # stat="probability",
)

sns.histplot(
    data=df_j,
    y="n_categories",
    hue="has_prophage",
    alpha=0.2,
    ax=axs[1, 1],
    legend=False,
    log_scale=True,
    bins=25,
    element="step",
    common_bins=True,
    # common_norm=False,
    # stat="probability",
)

axs[0, 1].set_visible(False)
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "joint_plot_prophages.png", dpi=300)
plt.show()

# %%
