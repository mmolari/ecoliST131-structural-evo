# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import numpy as np
import json

import pypangraph as pp
import utils as ut

from Bio import Phylo


def count_events(path_cats):
    """
    Takes as input paths that:
    - have two categories
    - the minority category is present in only one isolate
    Decides whether the minority category was subject to a gain or a loss.
    """
    assert len(path_cats) == 2
    p1, p2 = path_cats
    c1, n1, i1 = p1
    c2, n2, i2 = p2
    assert c2 == 1, f"{c1}, {c2}"
    assert n1[0] == n2[0], f"{n1}, {n2}"
    assert n1[-1] == n2[-1], f"{n1}, {n2}"
    n1, n2 = set(n1), set(n2)
    if n1 == n2:
        print("equal", n1, n2)
        return i2, "other", set()
    if n1.issubset(n2):
        return i2, "gain", n2 - n1
    elif n2.issubset(n1):
        return i2, "loss", n1 - n2
    else:
        print("other", n1, n2)
        return i2, "other", (n1 | n2) - (n1 & n2)


dset = "ST131_ABC"
fig_fld = pathlib.Path(f"figs/n1/{dset}")
fig_fld.mkdir(exist_ok=True, parents=True)


def svfig(name):
    # plt.savefig(fig_fld / f"{name}.svg")
    plt.savefig(fig_fld / f"{name}.png", facecolor="white")


fld = pathlib.Path(f"../../results/{dset}")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

# %%
df = pd.read_csv(df_file, index_col=0)

mask = df["n_iso"] > 50
mask &= df["n_categories"] == 2
df = df[mask]

sdf = df[df["singleton"]].copy()
Js = sdf.index.to_list()

# %%

tree = Phylo.read(tree_file, "newick")
tree.root_at_midpoint()
tree.ladderize()
terminal_len = {b.name: b.branch_length for b in tree.get_terminals()}

# %%
edf = pd.DataFrame.from_dict(terminal_len, orient="index", columns=["branch_length"])
edf["gain"] = 0
edf["loss"] = 0
edf["other"] = 0
sdf["event_type"] = None
gainloss_info = {}
for j in Js:
    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    paths = ut.pangraph_to_path_dict(pan)
    path_cat = ut.path_categories(paths)

    iso, tp, bs = count_events(path_cat)
    edf.loc[iso, tp] += 1
    sdf.loc[j, "event_type"] = tp
    gainloss_info[j] = {"iso": iso[0], "type": tp, "event_blocks": [b.id for b in bs]}

with open(f"data/{dset}/singleton_gainloss_info.json", "w") as f:
    json.dump(gainloss_info, f)

# %%
edf["ev_count"] = edf["gain"] + edf["loss"] + edf["other"]
edf["#events > 0"] = edf["ev_count"] > 0
edf.to_csv(f"data/{dset}/singleton_event_df.csv")
edf

# %%
# event rates
tot_terminal_branch_len = edf["branch_length"].sum()
tot_n_gains = edf["gain"].sum()
tot_n_losses = edf["loss"].sum()
tot_n_others = edf["other"].sum()
tot_n_events = tot_n_gains + tot_n_losses + tot_n_others

gain_rate = tot_n_gains / tot_terminal_branch_len
loss_rate = tot_n_losses / tot_terminal_branch_len
other_rate = tot_n_others / tot_terminal_branch_len
event_rate = tot_n_events / tot_terminal_branch_len

print(f"total terminal branch length: {tot_terminal_branch_len}")
print(f"total n. of gains: {tot_n_gains}")
print(f"total n. of losses: {tot_n_losses}")
print(f"total n. of other events: {tot_n_others}")
print(f"total n. of events: {tot_n_events}")
print(f"gain rate: {gain_rate}")
print(f"loss rate: {loss_rate}")
print(f"other rate: {other_rate}")
print(f"event rate: {event_rate}")


# %%

fig, ax = plt.subplots(1, 1, figsize=(3, 5))
et = sdf["event_type"].value_counts().to_dict()

colors = {
    "gain": "#51c46f",
    "loss": "#ed786b",
    "other": "#b8b8b8",
}

tot = sum(et.values())
bottom = 0
for k, v in et.items():
    p = ax.bar("signleton events", v, 0.8, bottom=bottom, label=k, color=colors[k])
    bottom += v
    ax.bar_label(p, label_type="center")

# put legend on the right side
ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)

sns.despine()
plt.tight_layout()
svfig("singleton_event_types")
plt.show()

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
fig, axs = plt.subplots(1, 2, figsize=(5, 8))

ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
ax.set_ylim(len(edf) + 1, 0)
ax.set_ylabel("isolates")

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
