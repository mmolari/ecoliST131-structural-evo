# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from Bio import Phylo
import json
import argparse
from collections import defaultdict

mb_df_file = "../../results/ST131_full/plasmids/mob/plasmid_summary.tsv"
tree_file = "../../results/ST131_full/pangraph/asm20-100-5-filtered-coretree.nwk"
pls_json = "../../config/datasets/ST131_full/plasmids.json"

tree = Phylo.read(tree_file, "newick")
tree.root_at_midpoint()
tree.ladderize()

with open(pls_json, "r") as f:
    pls_dict = json.load(f)

# %%
df = pd.read_csv(mb_df_file, sep="\t", index_col=0)

drop_cols = [
    "num_contigs",
    "size",
    # "gc",
    # "md5",
    # "rep_type(s)",
    "rep_type_accession(s)",
    # "relaxase_type(s)",
    "relaxase_type_accession(s)",
    # "mpf_type",
    "mpf_type_accession(s)",
    # "orit_type(s)",
    "orit_accession(s)",
    # "predicted_mobility",
    # "mash_nearest_neighbor",
    # "mash_neighbor_distance",
    # "mash_neighbor_identification",
    # "primary_cluster_id",
    # "secondary_cluster_id",
    "predicted_host_range_overall_rank",
    "predicted_host_range_overall_name",
    "observed_host_range_ncbi_rank",
    "observed_host_range_ncbi_name",
    "reported_host_range_lit_rank",
    "reported_host_range_lit_name",
    "associated_pmid(s)",
]
df = df.drop(columns=drop_cols)
df
# %%

relevant = [
    "rep_type(s)",
    "relaxase_type(s)",
    "orit_type(s)",
    "predicted_mobility",
]

for r in relevant:
    print(df[r].value_counts())

# %%

df["rep_type(s)"][df["relaxase_type(s)"].str.contains("MOBF")].str.contains(
    "IncF"
).value_counts()
# %%
df["rep_type(s)"][df["orit_type(s)"].str.contains("MOBF")].str.contains(
    "IncF"
).value_counts()


# %%
Frel = df["relaxase_type(s)"].str.contains("MOBF")
Fori = df["orit_type(s)"].str.contains("MOBF")

df["rep_type(s)"][Frel & Fori].str.contains("IncF").value_counts()

# %%

df[relevant][Frel & Fori].value_counts()

# %%
for r in relevant:
    X = sum([v.split(",") for v in df[r].values], [])
    print(np.unique(X, return_counts=True))
# %%

# colormap for predicted mobility
pm_color = {
    "conjugative": "red",
    "mobilizable": "sandybrown",
    "non-mobilizable": "white",
}

cmap = sns.color_palette("hls", 6)

# assign relaxase colors
rl_col = defaultdict(lambda: "lightgray")
rl_col |= {
    "MOBF": cmap[0],
    "-": "white",
    "MOBF,MOBP": cmap[1],
    "MOBP": cmap[2],
    "MOBQ": cmap[3],
    "MOBP,MOBP": cmap[4],
    "MOBF,MOBF": cmap[5],
}


# assign orit colors by abundance
# or_ct = df["orit_type(s)"].value_counts()
or_col = defaultdict(lambda: "lightgray")
or_col |= {
    "MOBF": cmap[0],
    "-": "white",
    "MOBP": cmap[1],
    "MOBQ,MOBQ": cmap[2],
    "MOBP,MOB_unknown": cmap[3],
    "MOBH": cmap[5],
}


I = list([l.name for l in tree.get_terminals()])

fig, axs = plt.subplots(
    1,
    2,
    figsize=(10, 1 + 0.08 * len(I)),
    sharey=True,
    gridspec_kw={"width_ratios": [1, 2]},
)

# plot tree
ax = axs[0]
Phylo.draw(tree, lambda x: "", do_show=False, axes=ax)

# plot plasmids MOB
ax = axs[1]
for i, iso in enumerate(I):
    y = i + 1
    x = 1
    if not iso in pls_dict:
        continue
    for pl in sorted(
        pls_dict[iso],
        key=lambda x: (
            df.loc[x]["predicted_mobility"],
            "IncF" not in df.loc[x]["rep_type(s)"],
            "MOBF" not in df.loc[x]["relaxase_type(s)"],
            "MOBF" not in df.loc[x]["orit_type(s)"],
        ),
    ):
        row = df.loc[pl].to_dict()
        # predicted mobility
        c = pm_color[row["predicted_mobility"]]
        ax.scatter(x, y, color=c, marker="s", edgecolor="gray")
        # predicted replication type
        if "IncF" in row["rep_type(s)"]:
            ax.scatter(x + 0.2, y, color="black", marker=".")
        # predicted relaxase
        c = rl_col[row["relaxase_type(s)"]]
        if c != "white":
            ax.scatter(x + 0.4, y, color=c, marker="x")
        # predicted oriT
        c = or_col[row["orit_type(s)"]]
        if c != "white":
            ax.scatter(x + 0.6, y, color=c, marker="+")

        x += 1

# generate legends
leg = []
for k, v in pm_color.items():
    leg.append(
        mpl.lines.Line2D(
            [], [], ls="", color=v, marker="s", markeredgecolor="black", label=k
        )
    )
leg.append(mpl.lines.Line2D([], [], ls="", color="black", marker=".", label="IncF"))
leg.append(mpl.lines.Line2D([], [], ls="", color="black", marker="x", label="relaxase"))
leg.append(mpl.lines.Line2D([], [], ls="", color="black", marker="+", label="oriT"))
leg_1 = ax.legend(handles=leg, loc="upper right")

leg = []
for k, v in rl_col.items():
    if v in ["lightgray", "white"]:
        continue
    leg.append(mpl.lines.Line2D([], [], ls="", color=v, marker="x", label=k))
leg_2 = ax.legend(handles=leg, loc="center right", title="relaxase")

leg = []
for k, v in or_col.items():
    if v in ["lightgray", "white"]:
        continue
    leg.append(mpl.lines.Line2D([], [], ls="", color=v, marker="+", label=k))
leg_3 = ax.legend(handles=leg, loc="lower right", title="oriT")
ax.add_artist(leg_1)
ax.add_artist(leg_2)
ax.set_xlabel("n. plasmids")

for ax in axs:
    for i in range(len(I)):
        ax.axhline(i + 1, color="lightgray", lw=0.5, zorder=-1)

# remove horizontal spacing between subplots
plt.tight_layout()
plt.subplots_adjust(wspace=0)

plt.savefig("figs/mob.png", dpi=300, facecolor="white")
plt.show()

# %%
