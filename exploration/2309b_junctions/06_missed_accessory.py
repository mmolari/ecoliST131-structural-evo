# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
import utils as ut
from collections import defaultdict


len_thr = 500
dset = "ST131_full"
edge_file = f"../../results/{dset}/backbone_joints/asm20-100-5/core_edges.csv"
pangraph_file = f"../../results/{dset}/pangraph/asm20-100-5-polished.json"

# %%
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
paths = ut.pangraph_to_path_dict(pan)
edges_df = pd.read_csv(edge_file)

# %%


def is_core(node_id):
    return (bdf.loc[node_id, "len"] >= len_thr) and bdf.loc[node_id, "core"]


# core-junctions dataframe
jdf = []
for iso, path in paths.items():
    junctions = ut.path_junction_split(path, is_core)
    for J in junctions:
        edge = J.flanking_edge()
        L = sum(bdf.loc[node.id, "len"] for node in J.center.nodes)
        jdf.append({"iso": iso, "edge": edge.to_unique_str_id(), "len": L})
jdf = pd.DataFrame(jdf)
jdf = jdf.pivot_table(index="iso", columns="edge", values="len")
jdf
# %%

# verify compatibility with edge file
old_edge_df = pd.read_csv(edge_file)
old_core_edges = old_edge_df[old_edge_df["count"] == len(paths)]["edge"].values
old_core_edges = set(
    [ut.Edge.from_str_id(e).to_unique_str_id() for e in old_core_edges]
)
new_core_edges = set(jdf.columns[jdf.notna().sum() == len(paths)])
assert old_core_edges == new_core_edges

# %%

# select core-edges
core_mask = jdf.notna().sum() == len(paths)
core_edges = jdf.columns[core_mask]
jdf_core = jdf[core_edges]

# select accessory-edges
acc_edges = jdf.columns[~core_mask]
jdf_acc = jdf[acc_edges]


# relative genome size
core_genome = jdf_core.sum().sum() / len(paths)
acc_genome = jdf_acc.sum().sum() / len(paths)
print(f"accessory genome in backbone: {core_genome/1e6:.2f} Mbp")
print(f"accessory genome out of backbone: {acc_genome/1e6:.2f} Mbp")
print(f"n. backbone edges {len(core_edges)}")
print(f"n. off-backbone edges {len(acc_edges)}")


# %%
def plot_core_len_vs_frequency(core_df):
    # evaluate average length of non-empty junctions
    avg_len = core_df[core_df > 0].mean(axis=0)
    avg_len.fillna(0, inplace=True)
    freq = (core_df > 0).sum(axis=0) / len(core_df.index)
    g = sns.jointplot(
        x=freq,
        y=avg_len,
        height=8,
        kind="hist",
        joint_kws={"bins": (50), "log_scale": (False, True)},
        marginal_kws={"bins": 50},
    )
    g.set_axis_labels("non-empty frequency", "average length of non-empty junctions")
    # add grid to the joint distribution
    g.ax_joint.grid(True, alpha=0.3)
    # set title
    g.fig.suptitle("backbone core edges survey")


plot_core_len_vs_frequency(jdf_core)
plt.tight_layout()
plt.savefig("figs/missed_accessory/backbone_edges_survey.png")
plt.show()

# %%


def plot_accessory_vs_core(df_acc, df_core):
    # evaluate average length of non-empty junctions
    avg_len_core = df_core.sum(axis=0) / df_core.shape[0]
    avg_len_core.fillna(0, inplace=True)
    avg_len_acc = df_acc.sum(axis=0) / df_acc.shape[0]
    avg_len_acc.fillna(0, inplace=True)
    plt.hist(
        [avg_len_core, avg_len_acc],
        bins=50,
        label=["backbone", "off-backbone"],
        log=True,
        histtype="stepfilled",
        alpha=0.3,
        edgecolor="k",
    )
    # plt.xscale("log")
    plt.legend()
    plt.xlabel("total amount of accessory genome per isolate (bp)")
    plt.ylabel("n. core-edges")


plot_accessory_vs_core(jdf_acc, jdf_core)
plt.tight_layout()
plt.savefig("figs/missed_accessory/core_edges_gap_lengths.png")
plt.show()

# %%


def plot_offbackbone_PA(dfa):
    consensus = jdf_acc.isna().mode().T[0]
    iso_mask = (consensus == jdf_acc.isna()).all(axis=1)
    sdf = dfa.loc[~iso_mask].notna()
    order = sdf.sum(axis=0).sort_values(ascending=False).index
    sdf = sdf[order]
    fig, ax = plt.subplots(figsize=np.array(sdf.shape) * 0.18 + np.array([3, 1]))
    sns.heatmap(sdf.T, cbar=True, ax=ax, cmap="GnBu", cbar_kws={"shrink": 0.5})
    # add title
    ax.set_title("presence/absence of off-backbone core-genome edges")


plot_offbackbone_PA(jdf_acc)
plt.tight_layout()
plt.savefig("figs/missed_accessory/offbackbone_PA.png")
plt.show()


# %%
def impute(dfa):
    idf = defaultdict(lambda: {"imputable": 0, "n_breakpoints": 0})
    full_strain_set = set(paths.keys())
    for col in dfa.columns:
        srs = dfa[col]
        strains_with = set(srs[srs.notna()].index)
        strains_without = full_strain_set - strains_with

        if len(strains_with) > len(strains_without):
            culprit = strains_without
        else:
            culprit = strains_with

        sequence = dfa[col].sum()
        for iso in culprit:
            idf[iso]["imputable"] += sequence / len(culprit)
            idf[iso]["n_breakpoints"] += 1
    idf = pd.DataFrame(idf).T
    return idf


def plot_impute_accessory(dfa):
    # impute missing accessory genome
    idf = impute(dfa)
    idf = idf.sort_values("imputable", ascending=False)

    fig, axs = plt.subplots(
        2,
        1,
        figsize=(idf.shape[0] * 0.2 + 1, 8),
        sharex=True,
    )
    ax = axs[0]
    h = idf["imputable"] / len(paths)
    bpl = ax.bar(range(len(idf)), h, color="C1")
    ax.bar_label(
        bpl,
        labels=[f"{np.round(x/1e3):.0f} kbp" for x in h],
        rotation=90,
        padding=5,
    )
    ax.set_ylabel("imputable accessory sequence loss (bp)")
    ax.set_yscale("log")

    ax = axs[1]
    bpl = ax.bar(range(len(idf)), idf["n_breakpoints"], color="C1")
    ax.bar_label(
        bpl,
        labels=[f"{x:.0f}" for x in idf["n_breakpoints"]],
        rotation=90,
        padding=5,
    )
    ax.set_ylabel("n. backbone breakpoints")
    ax.set_yscale("log")

    for ax in axs:
        ax.set_xticks(range(len(idf)))
        ax.set_xticklabels(
            idf.sort_values("imputable", ascending=False).index, rotation=90
        )
        ax.grid(alpha=0.3)

    fig.suptitle("off-backbone edges")
    sns.despine(fig)
    plt.tight_layout()
    plt.savefig("figs/missed_accessory/offbackbone_imputation.png")
    plt.show()


plot_impute_accessory(jdf_acc)

# %%


def plot_offbackbone_frequency(dfa):
    # bimodal frequency, rare and common
    notna = dfa.notna().sum() / len(dfa.index)
    fig, ax = plt.subplots(1, 1, figsize=(10, 3))
    bins = (np.arange(len(dfa.index) + 2) - 0.5) / len(dfa.index)
    plt.hist(notna, bins=bins, color="C1")
    # add minor xticks every 0.05
    ax.set_xticks(np.arange(0, 1.01, 0.05), minor=True)
    ax.set_xticks(np.arange(0, 1.01, 0.1), minor=False)

    plt.xlabel("frequency of off-backbone edges")
    plt.ylabel("n. edges")
    plt.tight_layout()
    plt.savefig("figs/missed_accessory/offbackbone_frequency.png")
    plt.show()


plot_offbackbone_frequency(jdf_acc)

# %%
