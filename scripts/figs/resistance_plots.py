import pathlib
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser(description="Plot metadata")
    parser.add_argument("--coregenome_tree", type=str)
    parser.add_argument("--card_df", type=str)
    parser.add_argument("--ncbi_df", type=str)
    parser.add_argument("--id_threshold", type=float)
    parser.add_argument("--outdir", type=str)
    return parser.parse_args()


def parse_element(x, thr):
    """
    Assign `.` = 0
    and `100.00;100;00;100.00` = 3
    """
    if x == ".":
        return 0
    x = str(x)
    els = x.split(";")
    ct = 0
    for el in els:
        if float(el) > thr * 100.0:
            ct += 1
    return ct


def load_res_df(fname, thr):
    df = pd.read_csv(fname, sep="\t", index_col=0)
    # transform filename to accession number
    df.index = df.index.map(lambda x: pathlib.Path(x).stem)
    df.drop(columns="NUM_FOUND", inplace=True)
    # parse elements
    df = df.applymap(lambda x: parse_element(x, thr))
    return df


def plot_resistance_vs_tree(tree, df):
    Y, X = df.shape

    iso_y = {l.name: n + 1 for n, l in enumerate(tree.get_terminals())}

    fig, axs = plt.subplots(
        1,
        2,
        figsize=(0.2 * X + 3, 1 + 0.08 * Y),
        sharey=True,
        gridspec_kw={"width_ratios": [1, X * 0.10]},
    )

    # plot tree
    ax = axs[0]
    Phylo.draw(tree, lambda x: "", do_show=False, axes=ax)

    # get colormap
    cmap = mpl.cm.get_cmap("rainbow")
    norm = mpl.colors.Normalize(vmin=1, vmax=df.max().max())
    mapp = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    # plot count matrix
    ax = axs[1]
    xticks = {}
    x = 0
    for col in df.columns:
        xticks[x] = col
        for iso in tree.get_terminals():
            y = iso_y[iso.name]
            ct = df.loc[iso.name, col]
            if ct > 0:
                c = cmap(norm(ct))
                ax.scatter(x, y, color=c, marker="s", edgecolor="black")
        ax.axvline(x, c="lightgray", zorder=-1, alpha=0.3)
        x += 1
    ax.set_xticks(list(xticks.keys()))
    ax.set_xticklabels(list(xticks.values()), rotation=90)
    ax.set_title("n. resistance genes")
    plt.colorbar(mapp, ax=ax, shrink=0.1)
    ax.grid(True, alpha=0.3, axis="x")

    for ax in axs:
        for n in range(Y):
            ax.axhline(n + 1, c="lightgray", zorder=-1, alpha=0.3)
        ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    return fig, axs


if __name__ == "__main__":
    args = parse_args()

    # create output path
    svfld = pathlib.Path(args.outdir)
    svfld.mkdir(exist_ok=True, parents=True)

    # load tree
    tree = Phylo.read(args.coregenome_tree, "newick")
    tree.root_at_midpoint()
    tree.ladderize()

    # for each resistance db
    df_fname = {
        "card": args.card_df,
        "ncbi": args.ncbi_df,
    }

    for df_type in ["card", "ncbi"]:
        # load resistance df
        df = load_res_df(df_fname[df_type], args.id_threshold)

        # plot
        fig, axs = plot_resistance_vs_tree(tree, df)
        plt.savefig(svfld / f"resistance_{df_type}.pdf")
