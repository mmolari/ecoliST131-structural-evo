import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--j_stats", type=str)
    parser.add_argument("--j_pangenome", type=str)
    parser.add_argument("--out_csv", type=str)
    parser.add_argument("--out_fig", type=str)
    parser.add_argument("--min_len", type=int)
    parser.add_argument("--min_paths", type=int)
    return parser.parse_args()


def load_data(j_stats, j_pangenome):
    df = pd.read_csv(j_stats, index_col=0)
    df2 = pd.read_csv(j_pangenome, index_col=0)
    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    mask = df["transitive"]
    df = df[~mask]
    return df


def make_figure(df, mask, min_len, min_paths, savename):
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    sns.scatterplot(
        data=df,
        x="pangenome_len",
        y="n_categories",
        hue=mask,
        legend=False,
        ax=ax,
        size=4,
    )
    plt.axvline(min_len, color="red", linestyle="--")
    plt.axhline(min_paths, color="red", linestyle="--")
    plt.xscale("log")
    plt.yscale("log")
    plt.text(0.1, 0.9, f"n. hotspots={mask.sum()}", transform=plt.gca().transAxes)
    plt.xlabel("pangenome length")
    plt.ylabel("number of paths")
    sns.despine()
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


if __name__ == "__main__":
    args = parse_args()

    df = load_data(args.j_stats, args.j_pangenome)

    hs_mask = df["pangenome_len"] >= args.min_len
    hs_mask &= df["n_categories"] >= args.min_paths

    make_figure(df, hs_mask, args.min_len, args.min_paths, args.out_fig)

    hs = df[hs_mask].sort_values(
        ["pangenome_len", "n_categories"],
        ascending=False,
    )
    hs.to_csv(args.out_csv)
