# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypangraph as pp
from Bio import Phylo
from collections import defaultdict
import pathlib

fig_fld = pathlib.Path("figs/f02")
fig_fld.mkdir(parents=True, exist_ok=True)

tree_fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_fname, "newick")
strains = [cl.name for cl in tree.get_terminals()]

fname = "../../results/ST131_ABC/hotspots/asm20-100-5/hotspots.csv"
hsdf = pd.read_csv(fname, index_col=0)

# hs = "XXVMWZCEKI_r__YUOECYBHUS_r"
# hs = "BWEZXGGFBK_r__MVMOFPVELT_r"
hs = "GPKQYOCEJI_r__NKVSUZGURN_f"

for hs in hsdf.index[::-1]:

    ML = hsdf.loc[hs, "max_length"] * 10 / 100000
    fig_size = (2 + max(10, ML), 12)
    width_ratios = [2, fig_size[0] - 2]
    print(ML)

    fname = (
        f"../../results/ST131_ABC/backbone_joints/asm20-100-5/joints_pangraph/{hs}.json"
    )
    pan = pp.Pangraph.load_json(fname)
    bdf = pan.to_blockstats_df()

    fig, axs = plt.subplots(
        # 1, 2, figsize=(12, 12), sharey=True, gridspec_kw={"width_ratios": [1, 4]}
        1,
        2,
        figsize=fig_size,
        sharey=True,
        gridspec_kw={"width_ratios": width_ratios},
    )

    # draw tree
    ax = axs[0]
    Phylo.draw(
        tree, axes=ax, label_func=lambda x: "", do_show=False, show_confidence=False
    )

    # draw pangraph
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

    block_colors = defaultdict(lambda: next(gen))
    ax = axs[1]
    for k, iso in enumerate(strains):
        y = k + 1
        if iso not in pan.strains():
            continue
        path = pan.paths[iso]
        Bs = path.block_ids
        Ss = path.block_strands
        x = 0
        for bid, strand in zip(Bs, Ss):
            color = block_colors[bid]
            l = bdf.loc[bid]["len"]
            edgecolor = "none" if strand else "k"
            ax.barh(
                y,
                l,
                left=x,
                height=0.8,
                color=color,
                edgecolor=edgecolor,
                linewidth=1,
            )
            x += l
    ax.set_xlabel("base pairs")
    pl = hsdf.loc[hs, "pangenome_len"]
    npaths = hsdf.loc[hs, "n_categories"]
    ax.set_title(f"{hs} | pangenome len = {pl//1000} kbp | n. paths = {npaths}")

    # despine
    for ax in axs:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(fig_fld / f"{hs}.png")
    plt.show()

# %%
