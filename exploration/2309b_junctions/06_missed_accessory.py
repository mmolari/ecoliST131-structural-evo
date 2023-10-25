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


def encompassing_edges(path, core_f):
    res = {}

    first_bid = None
    done = False
    L = len(path)
    for i, node in enumerate(path.nodes):
        if core_f(node.id):
            first_bid = node.id
            break
    edge = ut.Edge(path.nodes[i], None)
    link = []
    i += 1
    while not done:
        node = path.nodes[i]
        is_core = core_f(node.id)
        if is_core:
            edge.right = node
            res[edge] = link
            link = []
            edge = ut.Edge(node, None)
            if node.id == first_bid:
                done = True
        else:
            link.append(node)
        i += 1
        i %= L
    return res


def core_f(block_id):
    keep = bdf.loc[block_id, "core"]
    keep &= bdf.loc[block_id, "len"] >= len_thr
    return keep


links = {}
for iso, path in paths.items():
    links[iso] = encompassing_edges(path, core_f)
# %%

full_df = []
for iso, link in links.items():
    for edge, nodes in link.items():
        res = {}
        L = sum([bdf.loc[n.id, "len"] for n in nodes], start=0)
        res["edge"] = edge.to_unique_str_id()
        res["iso"] = iso
        res["len"] = L
        res["n_nodes"] = len(nodes)
        full_df.append(res)
full_df = pd.DataFrame(full_df)
df = full_df.pivot(index="iso", columns="edge", values="len")
col_ord = df.notna().sum().sort_values(ascending=False).index
df = df[col_ord]
df

# %%

core_edges_1 = df.columns[df.notna().sum() == len(paths)].to_numpy()
core_edges_2 = edges_df[edges_df["count"] == len(paths)]["edge"].to_numpy()

ce1 = set([ut.Edge.from_str_id(e) for e in core_edges_1])
ce2 = set([ut.Edge.from_str_id(e) for e in core_edges_2])

assert ce1 == ce2
# %%

df_core = df[core_edges_1].copy()
df_acc = df.drop(core_edges_1, axis=1).copy()
# %%
fig, ax = plt.subplots(figsize=(20, 30))
sns.heatmap(df_acc.notna(), cmap="Blues", ax=ax)
plt.show()


fig, ax = plt.subplots(figsize=(20, 30))
sns.heatmap(df_acc, cmap="rainbow", ax=ax)
plt.show()
# %%

sns.histplot(df_core.mean(), log_scale=True, bins=50)
# plt.xscale("log")
plt.show()
# %%
sns.histplot(x=(df_core == 0).sum(), y=df_core.mean(), log_scale=(False, True), bins=50)
# %%

n_zero_core = (df_core == 0).sum()
avg_nonzero_len_core = df_core[df_core > 0].mean()
avg_nonzero_len_core.fillna(0, inplace=True)

x_bins = np.linspace(0, len(paths) + 1, 50)
y_bins = np.logspace(
    np.log10(avg_nonzero_len_core.min()), np.log10(avg_nonzero_len_core.max()), 50
)
sns.jointplot(
    x=n_zero_core,
    y=avg_nonzero_len_core,
    kind="hist",
    joint_kws={
        "bins": (x_bins, y_bins),
        # "log_scale": (False, True),
    },
)
plt.yscale("log")
plt.show()
# %%
is_core = df.notna().sum() == len(paths)
is_core.value_counts()
# %%
core_acc_sequence = df.sum()[is_core].sum() / len(paths)
acc_acc_sequence = df.sum()[~is_core].sum() / len(paths)

print(f"Core-acc sequence: {core_acc_sequence/1e6} Mbp")
print(f"Acc-acc sequence: {acc_acc_sequence/1e6} Mbp")

# %%
from collections import defaultdict

imputable = defaultdict(float)
n_breakpoints = defaultdict(int)
full_strain_set = set(paths.keys())
for col in df.columns:
    srs = df[col]
    strains_with = set(srs[srs.notna()].index)
    strains_without = full_strain_set - strains_with

    if len(strains_with) > len(strains_without):
        culprit = strains_without
    else:
        culprit = strains_with

    sequence = df[col].sum()
    for iso in culprit:
        imputable[iso] += sequence / len(culprit)
        n_breakpoints[iso] += 1
imputable = pd.Series(imputable)
n_breakpoints = pd.Series(n_breakpoints)
imputable
# %%
fig, axs = plt.subplots(2, 1, figsize=(5, 10), sharex=True)

order = imputable.sort_values(ascending=False).index
ax = axs[0]
ax.plot(imputable[order] / len(paths), "o")
ax.set_yscale("log")

ax = axs[1]
ax.plot(n_breakpoints[order], "o")
ax.set_yscale("log")

for ax in axs:
    ax.set_xticks(range(len(imputable)))
    ax.set_xticklabels(imputable.sort_values(ascending=False).index, rotation=90)
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# %%
