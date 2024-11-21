# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gainloss_utils as ut
import pathlib
import json

fig_fld = pathlib.Path("figs/f05")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-alignment/filtered_corealignment_info_size.json"
with open(fname, "r") as f:
    aln_info = json.load(f)
aln_L = aln_info["polished aln size"]


def svfig(svname):
    pass
    plt.savefig(fig_fld / f"{svname}.png", dpi=200)
    plt.savefig(fig_fld / f"{svname}.svg")


# filenames
load_fld = pathlib.Path("../../results/ST131_ABC/rates/asm20-100-5")
fnames = {
    "cdf": {
        "terminal": load_fld / "terminal_coldspot.csv",
        "internal": load_fld / "nonsingleton_junct_df.csv",
    },
    "bdf": {
        "terminal": load_fld / "terminal_branches.csv",
        "internal": load_fld / "nonsingleton_branch_df.csv",
    },
    "coldspots": load_fld / "coldspot_df.csv",
    "events": load_fld / "merged_events_df.csv",
}

# %% terminal
cdf = ut.load_cdf(fnames["cdf"]["terminal"])
bdf = pd.read_csv(fnames["bdf"]["terminal"], index_col=0)
bdf["n_events"] = bdf["gain"] + bdf["loss"] + bdf["other"]
bdf["n_muts"] = bdf["branch_length"] * aln_L


def barplot_events(cdf, ax):
    vc = cdf[["event_type", "cat"]].value_counts()

    xlab = ["gain", "loss", "other"]

    for x, k in enumerate(xlab):
        y = 0
        for c in cdf["cat"].unique().sort_values():
            dy = vc.get((k, c), 0)
            kwargs = {"bottom": y, "color": ut.cat_colors[c]}
            if x == 0:
                kwargs["label"] = c
            ax.bar(x, dy, **kwargs)
            y += dy
        ax.text(x, y, y, ha="center", va="bottom")
    ax.legend()

    ax.set_xticks(range(len(xlab)))
    ax.set_xticklabels(xlab)
    ax.set_ylabel("n. events")


fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.5))
barplot_events(cdf, ax)
sns.despine()
plt.tight_layout()
# ax.set_title("terminal coldspots")
svfig("terminal_coldspots_count")
plt.show()


# %%

fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

ax = axs[0]
sns.histplot(
    data=bdf,
    x="branch_length",
    y="n_events",
    discrete=(False, True),
    ax=ax,
)
ax.set_xlabel("branch length")
ax.set_ylabel("n. events")

ax = axs[1]
g = sns.histplot(
    data=bdf,
    x="branch_length",
    hue=bdf["n_events"] > 0,
    element="step",
    ax=ax,
)
ax.set_xlabel("branch length")
ax.set_ylabel("n. branches")
# add legend
g.legend_.set_title("")
new_labels = ["no events", r"$\geq 1$" + " events"]
for t, l in zip(g.legend_.texts, new_labels):
    t.set_text(l)
sns.despine()
plt.tight_layout()
svfig("suppl_terminal_branches_vs_events")
plt.show()

# %%


fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

ax = axs[0]
sns.histplot(
    data=bdf,
    x="n_muts",
    y="n_events",
    discrete=(False, True),
    ax=ax,
)
ax.set_xlabel("n. substitutions")
ax.set_ylabel("n. events")

ax = axs[1]
g = sns.histplot(
    data=bdf,
    x="n_muts",
    hue=bdf["n_events"] > 0,
    element="step",
    ax=ax,
)
ax.set_xlabel("n. substitutions")
ax.set_ylabel("n. branches")
# add legend
g.legend_.set_title("")
new_labels = ["no events", r"$\geq 1$" + " events"]
for t, l in zip(g.legend_.texts, new_labels):
    t.set_text(l)
sns.despine()
plt.tight_layout()
svfig("suppl_terminal_muts_vs_events")
plt.show()

# %%
