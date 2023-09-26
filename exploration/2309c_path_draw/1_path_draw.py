# %%
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pypangraph as pp
import utils as ut
from Bio import Phylo

fld = pathlib.Path("../../results/ST131")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")

df = pd.read_csv(df_file, index_col=0)
df
# %%
for J in df.index:
    print(f"processing {J}")

    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{J}.json"
    pan = pp.Pangraph.load_json(pan_file)

    # preprocess
    bdf = pan.to_blockstats_df()

    paths = ut.pangraph_to_path_dict(pan)

    left = np.array([p.nodes[0] for p in paths.values()])
    right = np.array([p.nodes[-1] for p in paths.values()])

    rm = []
    for N in [left, right]:
        if np.all(N == N[0]):
            nid = N[0].id
            # print("All paths with the same node")
            if bdf.loc[nid]["core"]:
                # print("The node is core")
                rm.append(nid)

    # filter paths
    paths = ut.filter_paths(paths, lambda x: x not in rm)
    bdf.drop(rm, inplace=True)

    # plot image
    fig, axs = plt.subplots(
        1, 2, figsize=(10, 12), sharey=True, gridspec_kw={"width_ratios": [1, 3.5]}
    )

    # plot tree
    ax = axs[0]
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda x: "",
    )
    iso_y = {leaf.name: y + 1 for y, leaf in enumerate(tree.get_terminals())}
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))
    ax.grid(which="both", axis="y", alpha=0.5)
    # isolate y coordinate

    # assign colors
    colors = mpl.cm.get_cmap("nipy_spectral")(np.linspace(0, 1, len(bdf)))
    color_dict = {nid: c for nid, c in zip(bdf.sort_values("len").index, colors)}

    # plot paths
    ax = axs[1]
    for iso, path in paths.items():
        x = 0
        for node in path.nodes:
            nid = node.id
            l = bdf.loc[nid]["len"]
            core = bdf.loc[nid]["core"]
            dupl = bdf.loc[nid]["duplicated"]
            c = color_dict[nid]
            s = node.strand
            ax.barh(
                iso_y[iso],
                l,
                left=x,
                color=c,
                align="center",
                height=0.5 + 0.3 * dupl,
                edgecolor="gray" if s else "red",
                # hatch="///" if core else "",
            )
            x += l
    ax.set_xlabel("position (bp)")
    ax.set_title(J)

    for ax in axs:
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.savefig(f"figs/paths_simple/{J}.png", facecolor="white")
    # plt.show()
    plt.close(fig)
    # break

# %%
