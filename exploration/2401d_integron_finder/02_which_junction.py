# %%
import pathlib
import json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import location_utils as lu
from Bio import SeqIO

fig_fld = pathlib.Path("figs/f02")
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


# load integron dataframe
fname = data_fld / "integron_regions.csv"
df_int = pd.read_csv(fname, index_col=0)

# load junction info dataframes
fld = pathlib.Path(f"../../results/{dset}")
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

for integron, row in df_int.iterrows():
    pb, pe = row["beg"], row["end"]
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
                "integron": integron,
                "junction": j,
                "iso": iso,
                "in_core": in_core,
                "in_acc": in_acc,
                "ib": pb,
                "ie": pe,
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
                "integron": integron,
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
pp_to_j.to_csv(data_fld / "integrons_to_joints.csv")
rand_pp_to_j = pd.DataFrame(rand_pp_to_j)
pp_to_j
# %%

df_int["assigned"] = pp_to_j["integron"].value_counts()
df_int["assigned"].fillna(0, inplace=True)
df_int["assigned_rand"] = rand_pp_to_j["integron"].value_counts()
df_int["assigned_rand"].fillna(0, inplace=True)

sns.histplot(
    data=df_int,
    x="assigned",
    discrete=True,
    element="step",
    label="geNomad annotations",
)
sns.histplot(
    data=df_int,
    x="assigned_rand",
    discrete=True,
    element="step",
    label="shuffled annotations",
    alpha=0.5,
    color="gray",
)
plt.legend()
plt.xlabel("number of joints")
plt.ylabel("number of integron annotations")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "joints_per_integron.png", dpi=300)
plt.show()

# %%
df_j["n_integrons"] = pp_to_j["junction"].value_counts()
df_j["n_integrons"].fillna(0, inplace=True)
df_j["has_integron"] = df_j["n_integrons"] > 0
df_j["integron_freq"] = (
    pp_to_j.groupby("junction").apply(lambda x: len(np.unique(x["iso"])))
    / df_j["n_iso"]
)
df_j["integron_freq"].fillna(0, inplace=True)

df_j["n_integrons_rand"] = rand_pp_to_j["junction"].value_counts()
df_j["n_integrons_rand"].fillna(0, inplace=True)
df_j["has_integron_rand"] = df_j["n_integrons_rand"] > 0
df_j["integron_freq_rand"] = (
    rand_pp_to_j.groupby("junction").apply(lambda x: len(np.unique(x["iso"])))
    / df_j["n_iso"]
)
df_j["integron_freq_rand"].fillna(0, inplace=True)


sns.histplot(
    data=df_j, x="n_integrons", discrete=True, element="step", label="geNomad joints"
)
sns.histplot(
    data=df_j,
    x="n_integrons_rand",
    discrete=True,
    element="step",
    label="shuffled joints",
    alpha=0.5,
    color="gray",
)
plt.xscale("symlog")
plt.legend()
plt.xlabel("number of integron annotations")
plt.ylabel("number of joints")
plt.yscale("log")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "integron_per_joint.png", dpi=300)
plt.show()


# %%
sns.scatterplot(
    data=df_j,
    x="nonempty_freq",
    y="integron_freq",
    hue="has_integron",
    alpha=0.3,
)
plt.plot([0, 1], [0, 1], color="gray", linestyle="--", zorder=-3, alpha=0.2)
plt.xlabel("nonempty frequency")
plt.ylabel("integron frequency")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "integron_vs_nonempty_freq.png", dpi=300)
plt.show()

# sns.scatterplot(
#     data=df_j,
#     x="nonempty_freq",
#     y="integron_freq",
#     alpha=0.3,
#     label="geNomad annotations",
# )
# sns.scatterplot(
#     data=df_j,
#     x="nonempty_freq",
#     y="integron_freq_rand",
#     alpha=0.3,
#     label="shuffled annotations",
# )
# plt.plot([0, 1], [0, 1], color="gray", linestyle="--", zorder=-3, alpha=0.2)
# plt.xlabel("nonempty frequency")
# plt.ylabel("integron frequency")
# plt.legend()
# sns.despine()
# plt.tight_layout()
# plt.savefig(fig_fld / "integron_vs_nonempty_freq_rand.png", dpi=300)
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
    alpha=0.2,
    ax=ax,
)
sns.scatterplot(
    data=df_j[df_j["has_integron"]],
    x="delta_len",
    y="n_categories",
    hue="integron_freq",
    palette="viridis",
    ax=ax,
)
ax.set_xlabel("delta length")
ax.set_ylabel("n categories")

sns.histplot(
    data=df_j,
    x="delta_len",
    alpha=0.2,
    ax=axs[0, 0],
    log_scale=True,
    bins=25,
    element="step",
)

sns.histplot(
    data=df_j,
    y="n_categories",
    alpha=0.2,
    ax=axs[1, 1],
    log_scale=True,
    bins=25,
    element="step",
)

axs[0, 1].set_visible(False)
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "joint_plot_integrons.png", dpi=300)
plt.show()

# %%
