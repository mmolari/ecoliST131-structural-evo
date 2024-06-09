# %%
import numpy as np
import pandas as pd
import path_utils as pu
import pypangraph as pp
from Bio import Phylo
import matplotlib as mpl
import matplotlib.pyplot as plt
import pathlib
import argparse
import json
import itertools as itt


def show():
    # plt.show()
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hs", type=str)
    return parser.parse_args()


def load_data(hs):
    fname = f"data/{hs}/joint_graph.json"
    pan = pp.Pangraph.load_json(fname)

    fname = "data/tree.nwk"
    tree = Phylo.read(fname, "newick")
    return pan, tree


def draw_tree(tree, path_cats, sv_fld):
    leaves_pos = {t.name: i + 1 for i, t in enumerate(tree.get_terminals())}
    root_distance = {t.name: tree.distance(t) for t in tree.get_terminals()}

    fig, ax = plt.subplots(figsize=(3.5, 12))

    Phylo.draw(
        tree,
        axes=ax,
        show_confidence=False,
        label_func=lambda x: "",
        do_show=False,
    )
    for n, p, I in path_cats:
        y = [leaves_pos[i] for i in I]
        x = [root_distance[i] for i in I]
        c = "red" if n == 1 else None
        m = "x" if n == 1 else "o"
        ax.scatter(x, y, zorder=3, c=c, marker=m, s=15)
    ax.set_xticks([0, 5e-5, 1e-4])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(sv_fld / "path_categories.svg")
    plt.savefig(sv_fld / "path_categories.png", dpi=200)
    show()


def distance_df(pan):
    bdf = pan.to_blockstats_df()
    Ls = bdf["len"]
    padf = pan.to_blockcount_df()
    strains = pan.strains()
    res = []
    for i, j in itt.combinations(strains, 2):
        block_diff = np.abs(padf.loc[i] - padf.loc[j])
        delta_l = np.sum(block_diff * Ls)
        delta_n = np.sum(block_diff)
        res.append((i, j, delta_l, delta_n))
        res.append((j, i, delta_l, delta_n))
    dist_df = pd.DataFrame(res, columns=["i", "j", "delta_l", "delta_n"])
    dist_df.set_index(["i", "j"], inplace=True)
    return dist_df


def distance_plot(tree, dist_df, svfld):
    for lab, fname, title in [
        ("delta_l", "block_diffs_len.png", "len. private sequence"),
        ("delta_n", "block_diffs_num.png", "block n. difference"),
    ]:
        terminals = [t.name for t in tree.get_terminals()]
        N = len(terminals)
        M = np.zeros((N, N))
        for ni in range(N):
            for nj in range(ni + 1, N):
                i = terminals[ni]
                j = terminals[nj]
                if (i, j) in dist_df.index:
                    M[ni, nj] = dist_df.loc[i, j][lab]
                    M[nj, ni] = dist_df.loc[j, i][lab]
                else:
                    M[ni, nj] = np.nan
                    M[nj, ni] = np.nan

        # log norm
        Max = np.nanmax(M)
        Min = np.nanmin(M)
        if lab == "delta_l":
            norm = mpl.colors.LogNorm(vmin=max(Min, 100), vmax=Max, clip=True)
        else:
            norm = mpl.colors.Normalize(vmin=Min, vmax=Max)
        plt.figure(figsize=(10, 8))
        plt.imshow(M, norm=norm, cmap="coolwarm", interpolation="none")
        plt.title(title)
        # small colorbar
        plt.colorbar(shrink=0.5)
        plt.tight_layout()
        plt.savefig(svfld / fname, dpi=200)
        show()


if __name__ == "__main__":

    args = parse_args()
    hs = args.hs

    # hs = "CIRMBUYJFK_f__CWCCKOQCWZ_r"

    res_fld = pathlib.Path(f"res/{hs}/")
    res_fld.mkdir(exist_ok=True, parents=True)

    fig_fld = pathlib.Path(f"figs/{hs}/")
    fig_fld.mkdir(exist_ok=True, parents=True)

    pan, tree = load_data(hs)
    paths = pu.pangraph_to_path_dict(pan)
    path_cats = pu.path_categories(paths)

    # save path categories in json file
    fname = res_fld / "path_categories.json"
    with open(fname, "w") as f:
        jd = [(n, I) for n, path, I in path_cats]
        json.dump(jd, f, indent=2)

    draw_tree(tree, path_cats, fig_fld)

    ddf = distance_df(pan)
    distance_plot(tree, ddf, fig_fld)

# %%
