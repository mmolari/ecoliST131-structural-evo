# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib

data_fld = pathlib.Path("data/excluded_events")
data_fld.mkdir(exist_ok=True, parents=True)


def assign_mge_category(df):
    df["cat"] = None
    keys = ["genomad", "integrons", "isescan", "defensefinder"]
    F, T = False, True
    for val, k in [
        ((F, F, T, F), "IS"),
        ((F, F, F, F), "none"),
        ((T, F, T, F), "prophage"),
        ((T, F, F, F), "prophage"),
        ((F, T, T, T), "integron"),
        ((F, F, T, T), "defense"),
        ((T, F, F, T), "prophage"),
    ]:
        mask = np.all((df[keys] > 0) == val, axis=1)
        df.loc[mask, "cat"] = k

    # check that no "NaN" is left
    assert df["cat"].isna().sum() == 0, "some categories are not assigned"

    # make ordered categorical variable
    df["cat"] = pd.Categorical(
        df["cat"],
        categories=["IS", "integron", "prophage", "defense", "none"],
        ordered=True,
    )

    return df


cat_colors = {
    "IS": "C0",
    "prophage": "C4",
    "integron": "C1",
    "none": "#b8b8b8",
    "defense": "C2",
}


def plot_events(cdf, fname):

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
                "color": cat_colors[cat],
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
    plt.savefig(fname, dpi=150)
    plt.show()


fig_fld_main = pathlib.Path("figs/f04")

for dset in ["ST131_ABC", "ST131_sub_BC", "ST131_sub_C", "ST131_sub_C2"]:
    fig_fld = fig_fld_main / dset
    fig_fld.mkdir(exist_ok=True, parents=True)

    fname = f"../../results/{dset}/rates/asm20-100-5/nonsingleton_junct_df.csv"
    cdf = pd.read_csv(fname, index_col=0)
    fname = f"../../results/{dset}/rates/asm20-100-5/nonsingleton_branch_df.csv"
    bdf = pd.read_csv(fname, index_col=0)
    cdf["n_events"] = cdf["gain"] + cdf["loss"]
    cdf = assign_mge_category(cdf)

    plot_events(cdf, fig_fld / "internal_gainloss.png")

    mask = cdf["n_events"] <= 2
    cdf_correct = cdf[mask].copy()
    plot_events(cdf_correct, fig_fld / "internal_gainloss_correct.png")

    cdf[~mask]["n_events"].to_csv(data_fld / f"{dset}_n_events.csv")

# %%
