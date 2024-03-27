# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gainloss_utils as ut
import pathlib

fig_fld = pathlib.Path("figs/f05")
fig_fld.mkdir(exist_ok=True, parents=True)


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


def plot_events(cdf, bdf, fname):

    fig, axs = plt.subplots(1, 2, figsize=(8, 4))

    ax = axs[0]
    sns.histplot(data=cdf, x="n_minority", y="n_events", ax=ax, discrete=True)
    # plot diagonal
    x = np.arange(0, cdf["n_minority"].max())
    ax.plot(x, x, "--", lw=1, c="gray")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    # ax.set_aspect("equal")
    ax.set_xlabel("n. isolates with minority pattern")
    ax.set_ylabel("n. events")

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
    svfig(fname)
    plt.show()


plot_events(cdf, bdf, "suppl_internal_gainloss")


# %%
def plot_branches(df, fname):
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))

    df["n_gainloss"] = df["n_gain"] + df["n_loss"] + df["n_other"]
    mask = ~df["terminal"]
    sdf = df[mask]

    ax = axs[0]
    sns.histplot(
        sdf,
        x="branch_length",
        y="n_gainloss",
        ax=ax,
        discrete=(False, True),
    )
    ax.set_xlabel("branch length")
    ax.set_ylabel("n. events")

    ax = axs[1]
    g = sns.histplot(
        data=sdf,
        x="branch_length",
        hue=sdf["n_gainloss"] > 0,
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
    svfig(fname)
    plt.show()


plot_branches(bdf, "suppl_internal_branches")

# %%

mask = cdf["n_events"] <= 10

bdf_correct = bdf.copy()
for j, row in cdf[~mask].iterrows():
    for gb in row["gain_branches"].split("|"):
        bdf_correct.loc[gb, "n_gain"] -= 1
    for lb in row["loss_branches"].split("|"):
        bdf_correct.loc[gb, "n_loss"] -= 1

plot_branches(bdf_correct, "suppl_internal_branches_correct")

cdf_correct = cdf[mask].copy()
plot_events(cdf_correct, bdf_correct, "suppl_internal_gainloss_correct")
# %%
