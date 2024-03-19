# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import utils as ut
import pathlib

fig_fld = pathlib.Path("figs/n2")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/rates/asm20-100-5/nonsingleton_junct_df.csv"
cdf = pd.read_csv(fname, index_col=0)
fname = "../../results/ST131_ABC/rates/asm20-100-5/nonsingleton_branch_df.csv"
bdf = pd.read_csv(fname, index_col=0)
cdf["n_events"] = cdf["gain"] + cdf["loss"]

cdf = ut.assign_mge_category(cdf)


# %%
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
    plt.savefig(fig_fld / fname)
    plt.show()


plot_events(cdf, bdf, "internal_gainloss.png")


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
    ax.set_ylabel("n. internal branches")
    # add legend
    g.legend_.set_title("")
    new_labels = ["no events", ">0 events"]
    for t, l in zip(g.legend_.texts, new_labels):
        t.set_text(l)

    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / fname)
    plt.show()


plot_branches(bdf, "internal_branches.png")

# %%

mask = cdf["n_events"] <= 10

bdf_correct = bdf.copy()
for j, row in cdf[~mask].iterrows():
    for gb in row["gain_branches"].split("|"):
        bdf_correct.loc[gb, "n_gain"] -= 1
    for lb in row["loss_branches"].split("|"):
        bdf_correct.loc[gb, "n_loss"] -= 1

plot_branches(bdf_correct, "internal_branches_correct.png")

cdf_correct = cdf[mask].copy()
plot_events(cdf_correct, bdf_correct, "internal_gainloss_correct.png")
# %%
