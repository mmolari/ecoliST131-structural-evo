import argparse
import json
import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Diagnostic figures for recombination filter"
    )
    parser.add_argument("--idxs", help="info file with indices", type=str)
    parser.add_argument("--size", help="infor file with sizes", type=str)
    parser.add_argument("--fig_full", help="full output figure", type=str)
    parser.add_argument("--fig_reduced", help="reduced output figure", type=str)
    args = parser.parse_args()
    return args


def parse_info(idx_file, size_file):
    with open(size_file, "r") as f:
        info = json.load(f)

    params = {
        "window": info["window size"],
        "threshold": info["nsnps max"],
        "L": info["core aln size"],
    }

    with open(idx_file, "r") as f:
        info = json.load(f)

    idxs = {"all": info["idxs_all"], "removed": info["idxs_removed"]}
    idxs["kept"] = list(set(idxs["all"]) - set(idxs["removed"]))

    return params, idxs


def __plot_window(ax, params, idxs):
    L, window, threshold = params["L"], params["window"], params["threshold"]
    Aidxs, kept_idxs = idxs["all"], idxs["kept"]
    window_bins = np.arange(0, L + window + 1, step=window)
    ct1, _ = np.histogram(Aidxs, bins=window_bins)
    ct2, _ = np.histogram(kept_idxs, bins=window_bins)
    ctbins = np.arange(np.max(np.hstack([ct1, ct2])) + 3) - 0.5

    ax.hist(ct1, bins=ctbins, alpha=0.4, label="pre-filter", color="C0")
    ax.hist(ct2, bins=ctbins, alpha=0.4, label="post-filter", color="C2")
    ax.set_xlim(left=0)
    ax.axvline(threshold, c="k", ls=":", label="threshold")
    ax.set_yscale("log")
    ax.set_xscale("symlog")
    ax.set_title(
        f"window size = {window} bp, threshold > {threshold} neighbouring SNPs"
    )
    ax.legend(loc="upper right")
    ax.set_xlabel("n. core-genome alignment polymorphic positions per window")
    ax.set_ylabel("n. positions")


def __plot_hist(ax, params, idxs):
    L = params["L"]
    Aidxs, remove_idxs, kept_idxs = idxs["all"], idxs["removed"], idxs["kept"]
    bins = np.arange(L + 10001, step=10000)
    ax.hist(Aidxs, bins=bins, label="pre-filter")
    ax.hist(remove_idxs, bins=bins, label="removed")
    # ax.hist(kept_idxs, bins=bins, label="post-filter", histtype="step")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel("core genome alignment")
    ax.set_ylabel("SNPs per 10kbp")
    ax.set_title(
        f"n. SNPs in core alignment before / after filtering = {len(Aidxs)} / {len(kept_idxs)}"
    )


def __plot_cumulative(ax, params, idxs, plot_removed):
    L = params["L"]
    Aidxs, remove_idxs, kept_idxs = idxs["all"], idxs["removed"], idxs["kept"]
    kwargs = {
        "bins": np.arange(L + 10001, step=10000),
        "cumulative": True,
        "density": True,
        "histtype": "step",
    }
    ax.hist(Aidxs, label="pre-filter", color="C0", **kwargs)
    if plot_removed:
        ax.hist(remove_idxs, label="removed", color="C1", **kwargs)
    ax.hist(kept_idxs, label="post-filter", color="C2", **kwargs)
    ax.legend(loc="upper left")
    ax.set_xlabel("core genome alignment")
    ax.set_ylabel("cumul. distr. of SNPs")


def diagnostic_plot(params, idxs, filename):
    """Diagnostic plot to check the effect of recombination filtering."""

    fig, axs = plt.subplots(3, 1, sharex=False, figsize=(10, 8))

    __plot_window(axs[0], params, idxs)
    __plot_hist(axs[1], params, idxs)
    __plot_cumulative(axs[2], params, idxs, plot_removed=True)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close(fig)


def diagnostic_plot_reduced(params, idxs, filename):
    """Diagnostic plot to check the effect of recombination filtering."""

    fig, axs = plt.subplots(2, 1, sharex=False, figsize=(10, 6))

    __plot_hist(axs[0], params, idxs)
    __plot_cumulative(axs[1], params, idxs, plot_removed=False)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close(fig)


if __name__ == "__main__":
    args = parse_args()
    params, idxs = parse_info(args.idxs, args.size)
    diagnostic_plot(params, idxs, args.fig_full)
    diagnostic_plot_reduced(params, idxs, args.fig_reduced)
