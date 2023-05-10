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

fig, axs = plt.subplot_mosaic(
    """
AB
AB
""",
    figsize=(7, 6),
)
for ax, t, lab in [
    (axs["A"], tree_old, "before"),
    (axs["B"], tree_new, "after"),
]:
    Phylo.draw(
        t,
        do_show=False,
        axes=ax,
        label_func=lambda x: "",
    )
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_title(f"tree {lab}")
    for k in ["top", "left", "right"]:
        ax.spines[k].set_visible(False)
plt.tight_layout()
plt.savefig(fig_fld / f"tree_comparison.png")
plt.show()

# %%
