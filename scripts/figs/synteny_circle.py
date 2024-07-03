# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pypangraph as pp
from Bio import SeqIO
import pandas as pd
import pathlib
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pan", type=str)
    parser.add_argument("--ref", type=str)
    parser.add_argument("--fa", type=str)
    parser.add_argument("--len_thr", type=int)
    parser.add_argument("--block_colors", type=str)
    parser.add_argument("--mergers", type=str)
    parser.add_argument("--fig", type=str)
    return parser.parse_args()


def read_fasta(fname):
    with open(fname, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            return record.seq


def gc_skew(seq, bin_width):
    n = len(seq)
    bins = int(n / bin_width)
    skew = np.zeros(bins)
    for i in range(bins):
        start = i * bin_width
        end = (i + 1) * bin_width
        skew[i] = (seq[start:end].count("G") - seq[start:end].count("C")) / bin_width
    return skew


def circle_plot(pan, fname, ref, seq, len_thr, bc, mg):

    bdf = pan.to_blockstats_df()
    Bs = pan.paths[ref].block_ids
    offset = pan.paths[ref].block_positions[0]
    L_ref = sum([bdf["len"][b] for b in Bs])
    # polar coordinates
    fig, ax = plt.subplots(
        figsize=(6, 6),
        subplot_kw={"projection": "polar"},
    )

    def keep_f(bid):
        return bdf.loc[bid, "core"] and (bdf.loc[bid, "len"] > len_thr)

    factor = 2 * np.pi / L_ref
    x = offset * factor
    for bid in Bs:
        if not keep_f(bid):
            c = "gray"
            lw = 1
            # continue
        elif bid in mg:
            meta_block = mg[bid]
            c = bc.loc[meta_block, "color"]
            lw = 3
        else:
            c = bc.loc[bid, "color"]
            lw = 3
        l = bdf.loc[bid, "len"]
        # plt.plot([x, x + l], [0, 0], color=c, lw=lw)
        l *= factor
        plt.barh(y=10, width=l, height=lw, left=x, color=c, edgecolor="none")
        x += l

    # add skew
    bw = 20000
    y_set = 7.5
    y_skew = gc_skew(seq, bw)
    y_skew = y_skew * 12 + y_set
    x_skew = np.arange(len(y_skew), dtype=float) * bw
    x_skew *= factor
    # periodic boundary
    y_skew = np.concatenate([y_skew, y_skew[:1]])
    x_skew = np.concatenate([x_skew, x_skew[:1]])
    ax.plot(x_skew, y_skew, color="gray", lw=1)
    # densify
    x_skew_d = np.linspace(0, 2 * np.pi, 1000)
    y_skew_d = np.interp(x_skew_d, x_skew, y_skew, period=2 * np.pi)
    ax.fill_between(
        x_skew_d, y_skew_d, y_set, color="C0", alpha=0.5, where=y_skew_d > y_set
    )
    ax.fill_between(
        x_skew_d, y_skew_d, y_set, color="C3", alpha=0.5, where=y_skew_d < y_set
    )

    # set y axis
    ax.set_yticks([y_set])
    ax.set_yticklabels([""])
    ax.set_ylim([0, 10])

    xt = np.arange(0, L_ref, int(1e6))
    ax.set_xticks(xt * factor)
    ax.set_xticklabels([f"{x/1e6:.0f} Mb" for x in xt])
    ax.set_title(f"Synteny blocks in {ref}")
    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    plt.close()


if __name__ == "__main__":

    args = parse_args()

    # load arguments
    pan = pp.Pangraph.load_json(args.pan)
    seq = read_fasta(args.fa)
    bc = pd.read_csv(args.block_colors, index_col=0)
    mg = pd.read_csv(args.mergers, index_col=0)
    mg = mg["0"].to_dict()

    circle_plot(pan, args.fig, args.ref, seq, args.len_thr, bc, mg)
