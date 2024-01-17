# %%
import pathlib
import json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

fig_fld = pathlib.Path("figs/f04")
fig_fld.mkdir(exist_ok=True, parents=True)

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

# load junction dataframes
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_j = pd.read_csv(df_file, index_col=0)
df_j["delta_len"] = df_j["max_length"] - df_j["min_length"]

df_edge_len = fld / "backbone_joints/asm20-100-5/edge_len.csv"
df_el = pd.read_csv(df_edge_len, index_col=0)

pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
with open(pos_file) as f:
    jp = json.load(f)


# tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

# def extend_df(df, df_el):
#     new_df = {}
#     # FG = full graph
#     # frequency of edges
#     new_df["edge_freq_FG"] = df_el.notna().sum(axis=0) / df_el.shape[0]
#     new_df["nonempty_freq_FG"] = (df_el > 0).sum(axis=0) / df_el.notna().sum(axis=0)
#     new_df["n_iso_FG"] = df_el.notna().sum(axis=0)
#     new_df["nonempty_acc_len_FG"] = df_el.sum(axis=0) / (df_el > 0).sum(axis=0)
#     new_df = pd.DataFrame(new_df)
#     new_df.index.name = "edge"
#     df = df.join(new_df, how="outer").sort_values("n_iso_FG", ascending=False)
#     return df


# %%

#       B                        E
#       |------------------------|
#            b           e
#            |-----------|                   subset
#     |-----------------------------|        superset
#                       |----------------|   start
#    |--------------|                        end
#                                     |----| no


def within_no_loop(B, E, b, e):
    assert B <= E, f"{B=} > {E=}"
    assert b <= e, f"{b=} > {e=}"
    if b >= B and e <= E:
        return "subset"
    elif b <= B and e >= E:
        return "superset"
    elif b >= B and b <= E:
        return "start"
    elif e >= B and e <= E:
        return "end"
    else:
        return "no"


#   0   E                        B    L
#   ----|                        |-----
#            b           e
#            |-----------|                   no
#   |-----|                                  start
#                               |-----|      end
#   |-|                                      subset


def within_loop(B, E, b, e):
    assert B > E, f"{B=} < {E=}"
    assert b <= e, f"{b=} > {e=}"
    b_in = (b >= B) or (b <= E)
    e_in = (e >= B) or (e <= E)
    if b_in and e_in:
        return "subset"
    elif b_in:
        return "start"
    elif e_in:
        return "end"
    else:
        return "no"


def within(B, E, b, e):
    if B <= E:
        return within_no_loop(B, E, b, e)
    else:
        return within_loop(B, E, b, e)


pp_pos = {j: defaultdict(list) for j in jp}
for pp, row in df_ph.iterrows():
    bp, ep = row["start"], row["end"]
    iso = row["iso"]
    for j, pos_d in jp.items():
        if not iso in pos_d:
            continue
        cb, ab, ae, ce, strand = pos_d[iso]
        in_core = within(cb, ce, bp, ep)
        if in_core == "no":
            continue
        in_acc = within(ab, ae, bp, ep)
        pp_pos[j][iso].append((pp, in_core, in_acc, bp, ep))

# %%
has_pp = {
    j: np.any([np.any([x[2] != "no" for x in p]) for iso, p in pp_d.items()])
    for j, pp_d in pp_pos.items()
}

# %%
df_j["has_prophage"] = has_pp

# %%
sns.scatterplot(
    data=df_j,
    x="mean_length",
    y="nonempty_freq",
    hue="has_prophage",
    alpha=0.2,
)
plt.xscale("log")
plt.yscale("log")
plt.show()

sns.scatterplot(
    data=df_j,
    x="delta_len",
    y="n_categories",
    hue="has_prophage",
    alpha=0.2,
)
plt.xscale("log")
plt.yscale("log")
plt.show()

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
