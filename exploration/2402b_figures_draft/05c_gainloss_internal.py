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
cdf = ut.load_cdf(fnames["cdf"]["internal"])
bdf = pd.read_csv(fnames["bdf"]["internal"], index_col=0)
bdf["n_events"] = bdf["n_gain"] + bdf["n_loss"] + bdf["n_other"]
bdf["n_muts"] = bdf["branch_length"] * aln_L


def plot_events(cdf, fname):

    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

    ax = axs[0]
    g = sns.histplot(data=cdf, x="n_minority", y="n_events", ax=ax, discrete=True)

    # plot diagonal
    # x = np.arange(0, cdf["n_minority"].max())
    x = np.arange(0, cdf["n_events"].max())
    ax.set_yticks(range(0, x.max() + 1 + 2, 2))
    ax.plot(x, x, "--", lw=1, c="gray", label="diagonal")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    # ax.set_aspect("equal")
    ax.set_xlabel("n. isolates with minority pattern")
    ax.set_ylabel("n. events")
    # ax.legend(loc="center right")

    ax = axs[1]
    xlabs = ["gain", "loss", "other"]
    for x, et in enumerate(xlabs):
        y = 0
        for cat in cdf["cat"].unique().sort_values():
            mask = cdf["cat"] == cat
            dy = cdf[mask][et].sum()
            kwargs = {
                "color": ut.cat_colors[cat],
                "width": 0.8,
            }
            if x == 0:
                kwargs["label"] = cat
            ax.bar(x, dy, bottom=y, **kwargs)
            y += dy
        ax.text(x, y, y, ha="center", va="bottom")
    ax.legend()
    ax.set_ylabel("n. events")
    ax.set_xticks(range(len(xlabs)))
    ax.set_xticklabels(xlabs)

    sns.despine()
    plt.tight_layout()

    # colorbar in inset
    cax = fig.add_axes([0.40, 0.45, 0.02, 0.35])
    plt.colorbar(
        g.collections[0],
        cax=cax,
        orientation="vertical",
        label="n. non-singleton junctions",
    )

    svfig(fname)
    plt.show()


plot_events(cdf, "suppl_internal_gainloss")


# %%
def plot_branches(df, fname, xlab, xkey):
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

    df["n_gainloss"] = df["n_gain"] + df["n_loss"] + df["n_other"]
    mask = ~df["terminal"]
    sdf = df[mask]

    ax = axs[0]
    sns.histplot(
        sdf,
        x=xkey,
        y="n_gainloss",
        ax=ax,
        discrete=(False, True),
    )
    ax.set_xlabel(xlab)
    ax.set_ylabel("n. events")

    ax = axs[1]
    g = sns.histplot(
        data=sdf,
        x=xkey,
        hue=sdf["n_gainloss"] > 0,
        element="step",
        ax=ax,
    )
    ax.set_xlabel(xlab)
    ax.set_ylabel("n. branches")
    # add legend
    g.legend_.set_title("")
    new_labels = ["no events", r"$\geq 1$" + " events"]
    for t, l in zip(g.legend_.texts, new_labels):
        t.set_text(l)

    sns.despine()
    plt.tight_layout()
    svfig(fname)
    plt.show()


plot_branches(bdf, "suppl_internal_branches", "branch length", "branch_length")
plot_branches(bdf, "suppl_internal_branches_muts", "n. mutations", "n_muts")

# %%

lab = "tRNA dupl. event"
ut.cat_colors[lab] = "#888888"

mask = cdf["n_events"] >= 10
cdf_correct = cdf.copy()
cdf_correct["cat"] = pd.Categorical(
    cdf_correct["cat"],
    categories=["IS", "integron", "prophage", "defense", "none", lab],
    ordered=True,
)
cdf_correct.loc[mask, "cat"] = lab

bdf_correct = bdf.copy()
for j, row in cdf_correct[mask].iterrows():
    for gb in row["gain_branches"].split("|"):
        bdf_correct.loc[gb, "n_gain"] -= 1
    for lb in row["loss_branches"].split("|"):
        bdf_correct.loc[gb, "n_loss"] -= 1

plot_branches(bdf_correct, "suppl_internal_branches_correct")

plot_events(cdf_correct, "suppl_internal_gainloss_correct")
# %%
