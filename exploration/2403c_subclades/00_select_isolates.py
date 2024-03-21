# %%
from Bio import Phylo
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import pathlib
import json

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(fname, "newick")

fig_fld = pathlib.Path("figs/f00")
fig_fld.mkdir(exist_ok=True, parents=True)
data_fld = pathlib.Path("data/subcl_config")
data_fld.mkdir(exist_ok=True, parents=True)


def draw_tree(tree, subclade, ax):
    leaves = [x.name for x in subclade.get_terminals()]
    Phylo.draw(
        tree,
        do_show=False,
        axes=ax,
        label_func=lambda x: "x" if x.name in leaves else "",
        label_colors=lambda x: "black" if x == "" else "red",
    )
    ax.set_xticks(np.arange(10) * 1e-5)


# %%

subclades = {
    "BC": tree.root[1],
    "C": tree.root[1][1][1][1][1][1],
    "C2": tree.root[1][1][1][1][1][1][3],
}

fig, axs = plt.subplots(1, 3, figsize=(10, 10), sharey=True)
for i, (name, subclade) in enumerate(subclades.items()):
    draw_tree(tree, subclade, axs[i])
    axs[i].set_title(name)
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "subclades.png")
plt.show()

# %%
src_fld = pathlib.Path("../../config/datasets/ST131_ABC")
for sc, stree in subclades.items():
    print(f"Subclade: {sc} -----------------")

    # new fld
    sv_fld = data_fld / f"ST131_sub_{sc}"
    sv_fld.mkdir(exist_ok=True, parents=True)

    # leaves
    leaves = [x.name for x in stree.get_terminals()]

    # plasmids.json
    with open(src_fld / "plasmids.json") as f:
        plasmids = json.load(f)
    pl_acc_list = sum(plasmids.values(), [])
    print(
        f"n. plasmid keys: {len(plasmids)}",
        f"n. plasmids: {len(pl_acc_list)}",
        f"total: {len(set(pl_acc_list))+ len(leaves)}",
        sep=" | ",
    )
    plasmids = {k: v for k, v in plasmids.items() if k in leaves}
    pl_acc_list = sum(plasmids.values(), [])
    print(
        f"n. plasmid keys: {len(plasmids)}",
        f"n. plasmids: {len(pl_acc_list)}",
        f"total: {len(set(pl_acc_list))+ len(leaves)}",
        sep=" | ",
    )
    with open(sv_fld / "plasmids.json", "w") as f:
        json.dump(plasmids, f)

    # acc_info.csv
    ai = pd.read_csv(src_fld / "acc_info.csv")
    mask = ai["id"].isin(leaves)
    mask |= ai["id"].isin(pl_acc_list)

    print(f"n. acc_info.csv: {ai.shape[0]}")
    ai = ai[mask]
    print(f"n. acc_info.csv: {ai.shape[0]}")
    ai.to_csv(sv_fld / "acc_info.csv", index=False)

    # chromosomes.txt
    with open(src_fld / "chromosomes.txt") as f:
        chrom = f.read().splitlines()
    print(f"n. chromosomes: {len(chrom)}")
    chrom = [x for x in chrom if x in leaves]
    print(f"n. chromosomes: {len(chrom)}")
    with open(sv_fld / "chromosomes.txt", "w") as f:
        f.write("\n".join(chrom))
    # excluded.txt - empty file
    with open(sv_fld / "excluded.txt", "w") as f:
        f.write("")
    # metadata.csv
    md = pd.read_csv(src_fld / "metadata.csv", index_col=0)
    print(f"n. metadata: {md.shape[0]}")
    md = md.loc[leaves]
    print(f"n. metadata: {md.shape[0]}")
    md.to_csv(sv_fld / "metadata.csv")


# %%
