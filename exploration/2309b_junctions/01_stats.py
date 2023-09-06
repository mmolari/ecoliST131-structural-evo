# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import numpy as np

import pypangraph as pp
import utils as ut

from collections import defaultdict
from Bio import Phylo


def path_categories(paths):
    """Returns a list of touples, one per non-empty path, with the following info:
    (count, path, [list of isolates])"""
    iso_list = defaultdict(list)
    n_paths = defaultdict(int)
    nodes = {}
    for iso, path in paths.items():
        if len(path.nodes) > 0:
            n_paths[path] += 1
            iso_list[path].append(iso)
            nodes[path] = path.nodes

    # sort by count
    path_cat = [(count, nodes[path], iso_list[path]) for path, count in n_paths.items()]
    path_cat.sort(key=lambda x: x[0], reverse=True)
    return path_cat


def count_events(path_cats):
    assert len(path_cats) == 2
    p1, p2 = path_cats
    c1, n1, i1 = p1
    c2, n2, i2 = p2
    assert c2 == 1, f"{c1}, {c2}"
    assert n1[0] == n2[0], f"{n1}, {n2}"
    assert n1[-1] == n2[-1], f"{n1}, {n2}"
    if set(n1).issubset(n2):
        return i2, "gain"
    elif set(n2).issubset(n1):
        return i2, "loss"
    else:
        print(n1, n2)
        return i2, "other"


fig_fld = pathlib.Path("figs")


def svfig(name):
    plt.savefig(fig_fld / f"{name}.svg")
    plt.savefig(fig_fld / f"{name}.png")


fld = pathlib.Path("/home/marco/ownCloud/neherlab/code/pangenome-evo/results/ST131")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

# %%
df = pd.read_csv(df_file, index_col=0)
dfs = df[df["singleton"]]
Js = dfs.index.to_list()

# %%

tree = Phylo.read(tree_file, "newick")
tree.ladderize()
terminal_len = {b.name: b.branch_length for b in tree.get_terminals()}

# %%

edf = pd.DataFrame.from_dict(terminal_len, orient="index", columns=["branch_length"])
edf["gain"] = 0
edf["loss"] = 0
edf["other"] = 0
for j in Js:
    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    paths = ut.pangraph_to_path_dict(pan)
    path_cat = path_categories(paths)

    iso, tp = count_events(path_cat)
    edf.loc[iso, tp] += 1

# %%
edf["ev_count"] = edf["gain"] + edf["loss"] + edf["other"]
edf["#events > 0"] = edf["ev_count"] > 0
edf

# %%

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
sns.scatterplot(data=edf, x="branch_length", y="ev_count", alpha=0.5, ax=ax)
sns.despine()
plt.xlabel("terminal branch length")
plt.ylabel("n. of events on terminal branch")
plt.tight_layout()
svfig("len_vs_events_scatter")
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
x_bins = np.linspace(0, edf["branch_length"].max() * 1.1, 20)
y_bins = np.arange(edf["ev_count"].max() + 2) - 0.5
sns.histplot(
    data=edf,
    x="branch_length",
    y="ev_count",
    bins=(x_bins, y_bins),
    cbar=True,
    cbar_kws={"label": "n. of isolates"},
    ax=ax,
)
sns.despine()
plt.xlabel("terminal branch length")
plt.ylabel("n. of events on terminal branch")
plt.tight_layout()
svfig("len_vs_events_hist")
plt.show()

# %%

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
bins = np.linspace(0, edf["branch_length"].max() * 1.1, 20)
mask = edf["ev_count"] > 0
ax.hist(edf["branch_length"][mask], bins=bins, label="n. events > 0", alpha=0.3)
ax.hist(edf["branch_length"][~mask], bins=bins, label="n. events = 0", alpha=0.3)
ax.set_xlabel("terminal branch length")
ax.set_ylabel("n. of isolates")
ax.legend()
sns.despine()
plt.tight_layout()
svfig("len_vs_bool_events_hist")
plt.show()


# %%

# barplot for event types
fig, axs = plt.subplots(1, 2, figsize=(6, 6))

ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
ax.set_ylim(len(edf) + 1, 0)

ax = axs[1]
ax.barh(
    y=edf.index,
    width=edf["gain"],
    color="green",
    label="gain",
    left=edf["loss"] + edf["other"],
)
ax.barh(
    y=edf.index,
    width=edf["loss"],
    color="red",
    label="loss",
    left=edf["other"],
)
ax.barh(y=edf.index, width=edf["other"], color="grey", label="other")
ax.set_ylim(len(edf), -1)
ax.set_yticks([])
ax.set_xlabel("n. of events on terminal branch")

plt.xticks(rotation=90)
plt.legend()
plt.tight_layout()
sns.despine()
svfig("tree_vs_events")
plt.show()


# %%


# %%
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return "{p:.0f}%  ({v:d})".format(p=pct, v=val)

    return my_autopct


gains = edf["gain"].sum()
losses = edf["loss"].sum()
others = edf["other"].sum()
values = [gains, losses, others]
labels = ["gains", "losses", "others"]

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
ax.pie(
    values,
    labels=labels,
    autopct=make_autopct(values),
    # startangle=90,
    colors=["lightgreen", "lightcoral", "silver"],
)
ax.set_title("n. events on terminal branches")
fig.set_facecolor("white")
plt.tight_layout()
svfig("pie_events")
plt.show()
# %%
