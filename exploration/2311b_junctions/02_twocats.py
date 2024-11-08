# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib

import numpy as np
import os
import json

import pypangraph as pp
import utils as ut
import mugration_utils as mu


from collections import defaultdict
from Bio import Phylo


dset = "ST131_ABC"
fig_fld = pathlib.Path(f"figs/n2/{dset}")
fig_fld.mkdir(exist_ok=True, parents=True)


def svfig(name):
    # plt.savefig(fig_fld / f"{name}.svg")
    plt.savefig(fig_fld / f"{name}.png", facecolor="white")


fld = pathlib.Path(f"../../results/{dset}")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

# %%
df = pd.read_csv(df_file, index_col=0)
mask = df["n_categories"] == 2
mask &= df["n_iso"] > 50
mask &= ~df["singleton"]
sdf = df[mask].copy()
Js = sdf.index.to_list()
# %%

tree = Phylo.read(tree_file, "newick")
tree.ladderize()
terminal_len = {b.name: b.branch_length for b in tree.get_terminals()}

# %%
two_cat_isolates = {}
PA_df = {}
for j in Js:
    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    paths = ut.pangraph_to_path_dict(pan)
    path_cat = ut.path_categories(paths)

    if not len(path_cat) == 2:
        continue

    n1, p1, i1 = path_cat[0]
    n2, p2, i2 = path_cat[1]

    if n2 == 1:
        continue

    two_cat_isolates[j] = (p1, p2, i1, i2)
    PA_df[j] = {i: "A" for i in i1} | {i: "B" for i in i2}
PA_df = pd.DataFrame.from_dict(PA_df, orient="columns")
PA_df.fillna(".", inplace=True)

# %%
data_fld = pathlib.Path(f"data/{dset}/mugration/")
data_fld.mkdir(exist_ok=True, parents=True)

states_file = data_fld / "PA_states.csv"
PA_df.to_csv(states_file)
inf_file = data_fld / "mugration_inference.json"
# %%

# perform mugration inference
mu.mugration_inference(tree_file, states_file, data_fld, inf_file)

# %%
# load inference results
with open(inf_file, "r") as f:
    inf = json.load(f)

# %%

tree = Phylo.read(tree_file, "newick")
branch_len = {b.name: b.branch_length for b in tree.get_terminals()}
branch_len |= {b.name: b.branch_length for b in tree.get_nonterminals()}
branch_df = pd.DataFrame.from_dict(
    branch_len, orient="index", columns=["branch_length"]
)
branch_df["n_events"] = 0
branch_df["terminal"] = False
for b in tree.get_terminals():
    branch_df.loc[b.name, "terminal"] = True
event_info = {}
for k in two_cat_isolates:
    p1, p2, i1, i2 = two_cat_isolates[k]
    p1, p2 = set(p1), set(p2)
    info = {}
    if p1.issubset(p2):
        info["minority"] = "+"
    elif p2.issubset(p1):
        info["minority"] = "-"
    else:
        info["minority"] = "?"

    events = inf[k]["pa_pattern"]["events"]
    info["n_events"] = len(events)
    info["n_isolates"] = len(i2)
    info["gain"] = 0
    info["loss"] = 0
    for n, kind in events:
        branch_df.loc[n, "n_events"] += 1
        if info["minority"] == "?":
            continue
        elif info["minority"] == "+":
            if kind == "A|B":
                info["loss"] += 1
            elif kind == "B|A":
                info["gain"] += 1
        elif info["minority"] == "-":
            if kind == "A|B":
                info["gain"] += 1
            elif kind == "B|A":
                info["loss"] += 1

    event_info[k] = info
branch_df["n. events"] = "> 0"
branch_df.loc[branch_df["n_events"] == 0, "n. events"] = "= 0"
idf = pd.DataFrame.from_dict(event_info, orient="index")

# save dataframe
df_file = data_fld / "internal_branch_events.csv"
idf.to_csv(df_file)

idf
# %%
# %%


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return "{p:.0f}%  ({v:d})".format(p=pct, v=val)

    return my_autopct


gains = idf["gain"].sum()
losses = idf["loss"].sum()
others = idf["n_events"].sum() - gains - losses
values = [gains, losses]
labels = ["gains", "losses"]


fig, axs = plt.subplots(2, 2, figsize=(10, 8))

ax = axs[0, 0]
ax.pie(
    values,
    labels=labels,
    autopct=make_autopct(values),
    # startangle=90,
    colors=["lightgreen", "lightcoral", "silver"],
)
ax.set_title("n. events")
fig.set_facecolor("white")

ax = axs[1, 0]
bins = np.arange(0, idf["n_isolates"].max() + 1) + 0.5
sns.histplot(
    idf,
    x="n_isolates",
    y="n_events",
    bins=(bins, bins),
    cbar=True,
    cbar_kws={"label": "Number of two-category paths"},
    ax=ax,
)
ax.plot([0, bins.max()], [0, bins.max()], color="black", linestyle="--")
ax.set_xlabel("number of isolates with minority path")
ax.set_ylabel("number of events")

ax = axs[0, 1]
sdf = branch_df[~branch_df["terminal"]]
bins = (
    # np.linspace(0, sdf["branch_length"].max() * 1.05, 15),
    np.logspace(-7, -4, 15),
    np.arange(0, sdf["n_events"].max() + 2) - 0.5,
)
sns.histplot(
    sdf,
    x="branch_length",
    y="n_events",
    cbar=True,
    bins=bins,
    cbar_kws={"label": "Number of internal branches"},
    ax=ax,
)
ax.set_xlabel("branch length")
ax.set_ylabel("number of events")
ax.set_xscale("log")

ax = axs[1, 1]
sns.histplot(
    sdf,
    x="branch_length",
    bins=bins[0],
    ax=ax,
    hue=sdf["n. events"],
    alpha=0.4,
    element="step",
)
ax.set_xscale("log")
ax.set_xlabel("branch length")
ax.set_ylabel("number of internal branches")

plt.tight_layout()
sns.despine()
plt.savefig(fig_fld / "two_cat_events.png", dpi=300)
plt.savefig(fig_fld / "two_cat_events.pdf")
plt.show()
# %%

aln_len = 2427416
sbdf = branch_df[~branch_df["terminal"]]
tot_internal_branch_len = sbdf["branch_length"].sum()
tot_n_internal_events = sbdf["n_events"].sum()
rate = tot_internal_branch_len / tot_n_internal_events
per_mut = rate * aln_len

print(f"total internal branch length: {tot_internal_branch_len:.2e}")
print(f"total number of internal events: {tot_n_internal_events}")
print(f"event rate: one event every {rate:.2e}")
print(f"events every n. mutations: {per_mut:.3}")

# %%
