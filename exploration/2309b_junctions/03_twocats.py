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


fig_fld = pathlib.Path("figs")


fld = pathlib.Path("/home/marco/ownCloud/neherlab/code/pangenome-evo/results/ST131")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

# %%
df = pd.read_csv(df_file, index_col=0)
Js = df.index.to_list()
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

# %%
sdf = df.loc[two_cat_isolates.keys()]
fgfld = fig_fld / "two_cat_isolates"
src_fld = "../../figs/ST131/backbone_joints/asm20-100-5/joints_linear_plot"
fgfld.mkdir(exist_ok=True)
for k in two_cat_isolates:
    # copy all figures in one place
    os.system(f"cp {src_fld}/{k}.png {fgfld}/")

# %%
states_file = "data/mugration/PA_states.csv"
PA_df.to_csv(states_file)
svfld = "data/mugration"
inf_file = "data/mugration/mugration_inference.json"
# %%

# perform mugration inference
mu.mugration_inference(tree_file, states_file, svfld, inf_file)

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
    info = {}
    if len(p1) == 2 and len(p2) == 3:
        info["minority"] = "+"
    elif len(p1) == 3 and len(p2) == 2:
        info["minority"] = "+"
    else:
        info["minority"] = "+"

    events = inf[k]["pa_pattern"]["events"]
    info["n_events"] = len(events)
    info["n_isolates"] = len(i2)
    info["gain"] = 0
    info["loss"] = 0
    for n, kind in events:
        branch_df.loc[n, "n_events"] += 1
        if info["minority"] == "?":
            continue
        a = kind == "B|A"
        b = info["minority"] == "+"

        if (a and b) or (not a and not b):
            info["gain"] += 1
        else:
            info["loss"] += 1

    event_info[k] = info
branch_df["n. events"] = "> 0"
branch_df.loc[branch_df["n_events"] == 0, "n. events"] = "= 0"
idf = pd.DataFrame.from_dict(event_info, orient="index")
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

tot_internal_branch_len = branch_df[~branch_df["terminal"]]["branch_length"].sum()
tot_n_internal_events = branch_df[~branch_df["terminal"]]["n_events"].sum()

print(f"total internal branch length: {tot_internal_branch_len:.2e}")
print(f"total number of internal events: {tot_n_internal_events}")
print(
    f"event rate: {tot_n_internal_events / tot_internal_branch_len:.2e} events per branch length"
)

# %%
