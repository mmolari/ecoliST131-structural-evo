# %%
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np

# %%

ker = "asm20-100-5"
dset = "ST131"
tree_1_file = f"../../results/{dset}/gubbins/results/gubbins_{dset}.final_tree.tre"
tree_2_file = f"../../results/{dset}/gubbins/pan_{ker}_results/pan_gubbins_{dset}__{ker}.final_tree.tre"
tree_ref_file = f"../../results/{dset}/pangraph/asm20-100-5-filtered-coretree.nwk"

# %%
tree_1 = Phylo.read(tree_1_file, "newick")
tree_1.root_at_midpoint()
tree_1.ladderize()

tree_2 = Phylo.read(tree_2_file, "newick")
tree_2.root_at_midpoint()
tree_2.ladderize()

tree_ref = Phylo.read(tree_ref_file, "newick")
tree_ref.root_at_midpoint()
tree_ref.ladderize()

# %%
fig, ax = plt.subplots(1, 3, figsize=(30, 10))
Phylo.draw(tree_1, axes=ax[0], do_show=False)
Phylo.draw(tree_2, axes=ax[1], do_show=False)
Phylo.draw(tree_ref, axes=ax[2], do_show=False)
plt.show()


# %%
import copy


def draw_tree(tree, ax):
    # T = copy.deepcopy(tree)
    T = tree
    # assign y positions:
    for i, n in enumerate(T.get_terminals()):
        n.ypos = i
    for n in T.get_nonterminals(order="postorder"):
        child_ypos = [child.ypos for child in n.clades]
        n.ypos = np.mean(child_ypos)
    # assign x positions based on branch length:
    T.root.xpos = 0
    for n in T.get_nonterminals(order="preorder"):
        for c in n.clades:
            c.xpos = n.xpos + c.branch_length
    # draw:
    for n in T.get_nonterminals(order="preorder"):
        x = n.xpos
        y = n.ypos
        ax.scatter(x, y, color="r", marker=".")
        for child in n.clades:
            x_child = child.xpos
            y_child = child.ypos
            ax.plot([x, x, x_child], [y, y_child, y_child], color="k")


fig, ax = plt.subplots(1, 3, figsize=(20, 10))
draw_tree(tree_1, ax[0])
draw_tree(tree_2, ax[1])
draw_tree(tree_ref, ax[2])
plt.show()

# %%
import os

os.system(
    f"""
    treeknit {tree_ref_file} {tree_1_file} -o data/{dset}/ska --auspice-view
    treeknit {tree_ref_file} {tree_2_file} -o data/{dset}/pan --auspice-view
    treeknit {tree_1_file} {tree_2_file} -o data/{dset}/gub --auspice-view
    """
)

# """
# python3 scripts/figs/homoplasies.py \
#     --tree results/ST131_ABC/pangraph/asm20-100-5-coretree.nwk \
#     --aln results/ST131_ABC/pangraph/asm20-100-5-alignment/corealignment.fa \
#     --aln_info results/ST131_ABC/pangraph/asm20-100-5-alignment/corealignment_info.json \
#     --filt_tree results/ST131_ABC/gubbins/results/gubbins_ST131_ABC.final_tree.tre \
#     --filt_aln results/ST131_ABC/pangraph/asm20-100-5-alignment/filtered_corealignment.fa \
#     --filt_aln_info results/ST131_ABC/pangraph/asm20-100-5-alignment/filtered_corealignment_info_size.json \
#     --hist_fig test/hist.pdf \
#     --tree_fig test/tree.pdf \
#     --homoplasies_fig test/homoplasies.pdf
# """
