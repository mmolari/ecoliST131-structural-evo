# %%

import pathlib

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

from Bio import Phylo


fig_p = pathlib.Path("fig")
fig_p.mkdir(exist_ok=True)


def svfig(svname):
    for k in ["pdf", "png"]:
        plt.savefig(fig_p / f"{svname}.{k}", dpi=300)


# %%

tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, format="newick")
tree.ladderize()

# %%
dist_file = "../../results/ST131/distances/summary-asm20-100-5.csv"
dist_df = pd.read_csv(dist_file, index_col=[0, 1])
T_dist = dist_df.reset_index().pivot(
    index="si", columns="sj", values="core_div_filtered"
)
E_dist = dist_df.reset_index().pivot(index="si", columns="sj", values="edge_PA")
# E_dist = dist_df.reset_index().pivot(index="si", columns="sj", values="edge_sharing")
# E_dist = dist_df.reset_index().pivot(index="si", columns="sj", values="edge_PA_reduced")
B_dist = dist_df.reset_index().pivot(index="si", columns="sj", values="block_PA")

# %%
leaves_order = [l.name for l in tree.get_terminals()]
B_dist = B_dist.reindex(leaves_order, axis=0).reindex(leaves_order, axis=1)
E_dist = E_dist.reindex(leaves_order, axis=0).reindex(leaves_order, axis=1)
T_dist = T_dist.reindex(leaves_order, axis=0).reindex(leaves_order, axis=1)

# %%

import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


def distance_to_newick(dist_df):
    # leaves
    leaf_names = list(dist_df.index)
    num_leaves = len(leaf_names)

    # distance matrix
    distance_matrix = np.array(dist_df)

    # Perform nearest-neighbor joining
    cond_dist_mat = squareform(distance_matrix)
    linkage = hierarchy.linkage(cond_dist_mat, method="single")
    tree = hierarchy.to_tree(linkage, rd=False)

    # Create a dictionary to map leaf ids to names
    leaf_map = {i: leaf_names[i] for i in range(num_leaves)}

    # Traverse the tree and generate Newick format
    def traverse(node):
        if node.is_leaf():
            return leaf_map[node.id]
        else:
            left = traverse(node.get_left())
            right = traverse(node.get_right())
            return f"({left}:{node.dist/2},{right}:{node.dist/2})"

    # display tree with matplotlib
    fig, ax = plt.subplots(figsize=(10, 10))
    hierarchy.dendrogram(linkage, labels=leaf_names, ax=ax)
    plt.show()

    newick = traverse(tree) + ";\n"
    return newick, tree


for D_df, name in zip([B_dist, E_dist, T_dist], ["block", "edge", "core"]):
    newick_tree, tree = distance_to_newick(D_df)
    with open(f"data/{name}_tree.nwk", "w") as f:
        f.write(newick_tree)


# %%

os.system(
    f"""
    treeknit data/core_tree.nwk data/block_tree.nwk -o data/tk_cb_1 --auspice-view
    treeknit {tree_file} data/block_tree.nwk -o data/tk_cb_2 --auspice-view
    treeknit data/edge_tree.nwk data/block_tree.nwk -o data/tk_be --auspice-view
    """
)

# %%

# load tree
tree = Phylo.read("data/block_tree.nwk", format="newick")
# tree = Phylo.read("data/tk_cb/block_tree_resolved.nwk", format="newick")

# display tree
fig, ax = plt.subplots(figsize=(10, 10))
Phylo.draw(tree, axes=ax)
plt.show()

# %%
