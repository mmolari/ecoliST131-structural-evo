# %%

import pathlib
import utils as ut
from Bio import Phylo


data_fld = pathlib.Path("data")
data_fld.mkdir(exist_ok=True)

fig_fld = pathlib.Path("figs")
fig_fld.mkdir(exist_ok=True)


tree = Phylo.read(ut.tree_file, "newick")
for iso, node in enumerate(tree.get_nonterminals()):
    node.name = f"node_{iso:02d}"
Phylo.write(tree, ut.named_nodes_tree_file, "newick", format_branch_length="%.5e")

# %%
