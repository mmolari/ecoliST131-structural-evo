# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
import utils as ut

len_thr = 500
dset = "ST131_full"
edge_file = f"../../results/{dset}/backbone_joints/asm20-100-5/core_edges.csv"
pangraph_file = f"../../results/{dset}/pangraph/asm20-100-5-polished.json"

# %%
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
paths = ut.pangraph_to_path_dict(pan)
del pan
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


plot_core_len_vs_frequency(jdf_core)
plt.tight_layout()
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
        label=["core", "accessory"],
        log=True,
        histtype="stepfilled",
        alpha=0.3,
        edgecolor="k",
    )
    # plt.xscale("log")
    plt.legend()
    plt.xlabel("total amount of accessory genome per isolate (bp)")
    plt.show()


plot_accessory_vs_core(jdf_acc, jdf_core)
plt.tight_layout()
plt.show()

# %%


def plot_accessory_frequency():
    # bimodal frequency, rare and common
    pass


# %%
def plot_impute_accessory():
    # number and genome length
    pass
