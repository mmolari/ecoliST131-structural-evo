# %%
import numpy as np
import pandas as pd
import pypangraph as pp
import matplotlib.pyplot as plt
import seaborn as sns
import stats_utils as ut
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pan", type=str)
    parser.add_argument("--pair_dist", type=str)
    parser.add_argument("--out_csv", type=str)
    parser.add_argument("--out_fig", type=str)
    return parser.parse_args()


def load_pair_dist(fname):
    ddf = pd.read_csv(fname)
    mask = ddf["si"] > ddf["sj"]
    ddf = ddf[mask]
    ddf.set_index(["si", "sj"], inplace=True)
    return ddf


def make_fig(df, fig_fname):

    fig, axs = plt.subplots(3, 3, figsize=(12, 12), sharex=True)

    for nax, k in enumerate(
        [
            "n. private blocks",
            "merge_n_private",
            "len. private blocks",
            "n. breakpoints",
            "local_core_div",
            "merge_n_edges",
            "merge_core_len",
            "merge_acc_len",
            "merge_n_shared",
        ]
    ):
        ax = axs.flatten()[nax]
        sns.scatterplot(
            data=df,
            x="core_div_filtered",
            y=k,
            alpha=0.1,
            marker=".",
            ax=ax,
        )
        sns.lineplot(
            data=df,
            x="core_div_filtered",
            y=k,
            alpha=0.5,
            color="k",
            estimator=np.mean,
            # errorbar="sd",
            ax=ax,
        )

        # group by x in bins of size 0.0001
        df["bin"] = pd.cut(df["core_div_filtered"], 11)
        # x coordinate
        xc = df.groupby("bin", observed=True)["core_div_filtered"].mean()
        yc = df.groupby("bin", observed=True)[k].mean()
        err = df.groupby("bin", observed=True)[k].std()
        ax.errorbar(xc, yc, yerr=err, fmt="o", color="C1", ls="--")
    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fname)
    plt.close(fig)


if __name__ == "__main__":
    args = parse_args()
    ddf = load_pair_dist(args.pair_dist)
    pan = pp.Pangraph.load_json(args.pan)
    df = ut.extract_hotspot_stats(pan)
    df["core_div_filtered"] = ddf.loc[df.index, "core_div_filtered"]
    df.to_csv(args.out_csv)
    make_fig(df, args.out_fig)
