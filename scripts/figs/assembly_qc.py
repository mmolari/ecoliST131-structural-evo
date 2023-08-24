import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="summary figure for BUSCO assembly quality control."
    )
    parser.add_argument("--csv", type=str, help="input summary csv")
    parser.add_argument("--fig", type=str, help="output figure")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # load dataframe
    df = pd.read_csv(args.csv, index_col="iso")
    N = len(df)

    fig, axs = plt.subplots(1, 2, figsize=(8, 15), sharey=True)

    # number of contigs
    ax = axs[0]
    y_iso = {iso: i for i, iso in enumerate(df.index[::-1])}
    y_ticks = [y_iso[iso] for iso in df.index]
    y_labels = [iso for iso in df.index]

    for iso in df.index:
        row = df.loc[iso]
        x = row["Number of contigs"]
        ax.barh(y_iso[iso], x, color="black")

    ax.set_xscale("log")
    ax.set_xlabel("number of contigs")
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_ylim(-1, N)
    ax.grid(alpha=0.3)

    # missing/fragmented
    ax = axs[1]

    colors = {
        # "Single copy": "green",
        # "Multi copy": "darkgreen",
        "Fragmented": "orange",
        "Missing": "red",
    }

    for iso in df.index:
        row = df.loc[iso]
        x = 0
        for k, c in colors.items():
            dx = row[k]
            ax.barh(
                y_iso[iso],
                dx,
                left=x,
                color=c,
                label=k if y_iso[iso] == 0 else None,
            )
            x += dx
    ax.set_xlabel("percent of BUSCO genes")
    ax.legend()
    ax.grid(alpha=0.3)
    sns.despine()
    plt.tight_layout()
    plt.savefig(args.fig, facecolor="white")
    plt.close(fig)
