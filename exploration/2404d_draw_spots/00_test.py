# %%
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import SeqIO, Phylo
import pypangraph as pp
from collections import defaultdict
import plot_utils as put


def interesting_coldspots():
    df = put.load_coldspots_df()
    mask = ~df["singleton"]
    df = df[mask].sort_values("majority_category")
    return df


cdf = interesting_coldspots()

coldspots = [
    "CYGJWOEQKN_f__SKPHAXSFLS_f",
    "OURQVJZAZZ_f__UTYAQKFQDH_f",
    "IHERCMJOSU_f__JJRRWBDVGH_f",  # good IS example
    "YVEVUPDYEE_f__ZFYFAGFPQE_f",  # good prophage example
    # "RYYAQMEJGY_r__ZTHKZYHPIX_f",
]

edge = coldspots[-1]
tree, strains = put.load_tree()
pan = put.load_pangraph(edge)
mges = put.load_MGEs(edge)


# %%


def plot_tree(ax, tree):
    Phylo.draw(
        tree,
        axes=ax,
        label_func=lambda x: "",
        do_show=False,
        show_confidence=False,
    )


def color_dict():

    def color_generator():
        cm = plt.get_cmap("tab20")
        for i in range(20):
            yield cm(i)
        cm = plt.get_cmap("tab20b")
        for i in range(20):
            yield cm(i)
        cm = plt.get_cmap("tab20c")
        for i in range(20):
            yield cm(i)
        cm = plt.get_cmap("hsv")
        while True:
            yield cm(np.random.rand())

    gen = color_generator()
    return defaultdict(lambda: next(gen))


def pan_block_lengths(pan):
    Ls = defaultdict(dict)
    for b in pan.blocks:
        bid = b.id
        aln = b.alignment
        S, O = aln.generate_sequences()
        for s, o in zip(S, O):
            Ls[bid][o] = len(s)
    return Ls


def plot_paths(ax, pan, strains):

    Bcolors = color_dict()
    # color first/last blocks:
    start_color = "#2db7a7ff"
    end_color = "#cf4b36ff"
    Bs = pan.paths[strains[0]].block_ids
    Bcolors[Bs[0]] = start_color
    Bcolors[Bs[-1]] = end_color

    bdf = pan.to_blockstats_df()
    Ls = pan_block_lengths(pan)
    skipped = set(strains) - set(pan.strains())

    for k, iso in enumerate(strains):
        if iso in skipped:
            continue
        path = pan.paths[iso]
        Bs = path.block_ids
        Os = path.block_nums
        Ss = path.block_strands

        y = k + 1
        x = 0
        for b, o, s in zip(Bs, Os, Ss):
            l = Ls[b][(iso, o, s)]
            c = Bcolors[b]
            ax.barh(y, l, height=0.8, left=x, color=c)
            x += l


def add_mges(ax, mges, strains):
    mge_color = {
        "IS": "red",
        "prophage": "yellow",
        "defense system": "green",
        "integron": "k",
    }

    y_iso = {iso: k + 1 for k, iso in enumerate(strains)}
    for idx, row in mges.iterrows():
        y = y_iso[row["iso"]]
        strand = row["j_strand"]
        if strand:
            start = row["jcb"]
            b = row["ib"] - start
            e = row["ie"] - start
        else:
            end = row["jce"]
            b = end - row["ie"]
            e = end - row["ib"]
        c = mge_color[row["MGE_type"]]
        ax.plot([b, e], [y, y], "|-", color=c)


fig, axs = plt.subplots(
    1, 2, sharey=True, figsize=(12, 12), gridspec_kw={"width_ratios": [1, 5]}
)

ax = axs[0]
plot_tree(ax, tree)
ax = axs[1]
plot_paths(ax, pan, strains)
add_mges(ax, mges, strains)

plt.tight_layout()
plt.show()

# %%


# %%
