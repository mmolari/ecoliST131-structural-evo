# %%
import pandas as pd
import matplotlib.pyplot as plt
import pathlib

from Bio import Phylo

fig_fld = pathlib.Path("figs")
fig_fld.mkdir(exist_ok=True)

tree_old = Phylo.read(
    "../../results/ST131/pangraph/asm20-100-5-coretree.nwk", format="newick"
)
tree_old.ladderize()

tree_new = Phylo.read(
    "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk", format="newick"
)
tree_new.ladderize()

for t,lab in [(tree_old, "before"), (tree_new, "after")]:
    fig, ax = plt.subplots(1, 1, figsize=(5,7))
    Phylo.draw(
        t,
        do_show=False,
        axes=ax,
        label_func=lambda x: "",
    )
    ax.set_ylabel("")
    ax.set_yticks([])
    for k in ["top", "left", "right"]:
        ax.spines[k].set_visible(False)
    plt.tight_layout()
    plt.savefig(fig_fld / f"tree_{lab}.png")
    plt.show()

# %%
