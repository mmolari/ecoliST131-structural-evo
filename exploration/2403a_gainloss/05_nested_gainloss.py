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


fname = "../../results/ST131_ABC/rates/asm20-100-5/AB_states.csv"
abs_df = pd.read_csv(fname, index_col=0)
fname = "../../results/ST131_ABC/rates/asm20-100-5/AB_nestedness.csv"
abn_df = pd.read_csv(fname, index_col=0)
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

    color_downstream(tree.root, "silver")

    for clade in tree.get_nonterminals() + tree.get_terminals():
        if clade.name in sdf["branch"].values:
            row = sdf[sdf["branch"] == clade.name].iloc[0]
            event = row["type"]
            color = colors[event]
            # print(f"color {clade.name}, {event}, {color}")
            color_downstream(clade, color)


def draw_PA(j, x, strains, ax):
    pa_j = abs_df[j]
    nst = abn_df.loc[j, "event_type"]
    match nst:
        case "A>B":
            pa_dict = {"A": True, "B": False}
        case "A<B":
            pa_dict = {"A": False, "B": True}
        case _:
            raise ValueError(f"Unknown nestedness {nst}")
    for n, iso in enumerate(strains):
        pa = pa_j[iso]
        present = pa_dict[pa]
        # filled or empty circle
        marker = "o" if present else "x"
        # color
        color = "k" if present else "lightgray"
        ax.scatter(x, n, color=color, marker=marker, s=20)


def draw_tree(sdf, ax):
    tree = load_tree()
    color_tree(tree, sdf)
    j = sdf["junction"].iloc[0]

    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: None)
    ax.set_title(j, size=10)

    strains = [x.name for x in tree.get_terminals()]
    x = max([tree.distance(x) for x in tree.get_terminals()])
    draw_PA(j, x, strains, ax)


for idx, sdf in edf[edf.junction.isin(two_events)].groupby("junction"):
    fig, ax = plt.subplots(1, 1, figsize=(4, 10), sharex=True, sharey=True)
    k = 0
    draw_tree(sdf, ax)
    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / f"{idx}.png", dpi=150)
    plt.show()
    # break

# %%
