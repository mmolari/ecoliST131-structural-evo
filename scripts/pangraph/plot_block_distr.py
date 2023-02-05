import pypangraph
import argparse
import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="plot block frequency and length cumulative distributions."
    )
    parser.add_argument("--pangraph", help="json pangraph file.", type=str)
    parser.add_argument("--fig", help="output figure", type=str)
    args = parser.parse_args()
    return args


def summary_plot(df, savename=None):

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    n = df["n. strains"]
    L = df["len"]
    core = df["core"]
    dupl = df["duplicated"]

    C1, C2 = "#0E7395", "#E11444"

    bins = np.arange(n.max() + 2) - 0.5
    kw = {"cumulative": True, "histtype": "step", "alpha": 0.7}

    # block frequency
    ax = axs[0]
    ax.hist(n, bins=bins, color=C1, **kw)
    ax.set_ylabel("n. blocks", color=C1)

    axt = ax.twinx()
    axt.hist(n, weights=L, bins=bins, color=C2, **kw)
    axt.set_ylabel("block length (bp)", color=C2)

    ax.set_xlabel("n. strains")
    ax.set_title("cumulative frequency distr.")
    ax.set_xlim(bins.min(), bins.max())

    # block length
    bins = np.logspace(1, np.log10(L.max()) + 0.2, 100)
    ax = axs[1]
    ax.hist(L, bins=bins, color=C1, **kw)
    ax.set_ylabel("n. blocks", color=C1)

    axt = ax.twinx()
    axt.hist(L, weights=L, bins=bins, color=C2, **kw)
    axt.set_ylabel("block length (bp)", color=C2)

    ax.set_title("cumulative length distr.")
    ax.set_xlabel("block length")
    ax.set_xscale("log")
    ax.set_xlim(bins.min(), bins.max())

    for ax in axs:
        ax.grid(alpha=0.3)

    # summary stats
    txts = [
        f"{len(n)} blocks -> {L.sum()/1e6:.3} Mbp",
        f"{core.sum()} core blocks -> {L[core].sum()/1e6:.3} Mbp",
        f"{dupl.sum()} dupl blocks -> {L[dupl].sum()/1e6:.3} Mbp",
    ]
    ax = axs[0]
    for idx, txt in enumerate(txts):
        ax.text(0.1, 0.9 - 0.05 * idx, txt, transform=ax.transAxes)

    # finalize and save/show
    plt.tight_layout()
    if savename is None:
        plt.show()
    else:
        plt.savefig(savename)
        plt.close(fig)


if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # load pangraph
    pan = pypangraph.Pangraph.load_json(args.pangraph)

    # produce stats dataframe
    df = pan.to_blockstats_df()

    # produce and save plot
    summary_plot(df, savename=args.fig)
