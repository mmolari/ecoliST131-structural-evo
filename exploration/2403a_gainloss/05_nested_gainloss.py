# %%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
from Bio import Phylo

fig_fld = pathlib.Path("figs/n5")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/rates/asm20-100-5/merged_events_df.csv"
edf = pd.read_csv(fname)
fname = "../../results/ST131_ABC/rates/asm20-100-5/nonsingleton_branch_df.csv"
bdf = pd.read_csv(fname, index_col=0)


L_tot = bdf["branch_length"].sum()
L_ter = bdf[bdf["terminal"]]["branch_length"].sum()
L_int = bdf[~bdf["terminal"]]["branch_length"].sum()

# excluded junctions
j_count = edf["junction"].value_counts()
excl_j = j_count[j_count > 10].index

# %%
non_singleton = edf[~edf.singleton]
print(non_singleton.junction.value_counts().value_counts())

# %%
one_event = non_singleton.junction.value_counts() == 1
one_event = one_event[one_event].index

edf[edf.junction.isin(one_event)]["type"].value_counts()

# %%
two_events = non_singleton.junction.value_counts() == 2
two_events = two_events[two_events].index

# %%


def load_tree():
    fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
    tree = Phylo.read(fname, "newick")

    # assign names to internal nodes
    for k, node in enumerate(tree.get_nonterminals()):
        node.name = f"int_node_{k}"

    return tree


def color_downstream(branch, color):
    branch.color = color
    for clade in branch.clades:
        color_downstream(clade, color)


def color_tree(tree, sdf):
    colors = {
        "gain": "green",
        "loss": "red",
        "other": "purple",
    }

    for clade in tree.get_nonterminals() + tree.get_terminals():
        if clade.name in sdf["branch"].values:
            row = sdf[sdf["branch"] == clade.name].iloc[0]
            event = row["type"]
            color = colors[event]
            # print(f"color {clade.name}, {event}, {color}")
            color_downstream(clade, color)


def draw_tree(sdf, ax):

    tree = load_tree()
    color_tree(tree, sdf)

    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: None)
    ax.set_title(sdf["junction"].iloc[0], size=10)


fig, axs = plt.subplots(1, 4, figsize=(10, 10), sharex=True, sharey=True)
k = 0
for idx, sdf in edf[edf.junction.isin(two_events)].groupby("junction"):

    if sdf["type"].nunique() == 1:
        continue

    ax = axs[k]
    draw_tree(sdf, ax)
    k += 1
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "two_events.png", dpi=150)
plt.show()

# %%
