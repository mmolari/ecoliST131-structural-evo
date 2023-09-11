# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import numpy as np

import pypangraph as pp
import utils as ut

from collections import defaultdict
from Bio import Phylo


fig_fld = pathlib.Path("figs")


def svfig(name):
    # plt.savefig(fig_fld / f"{name}.svg")
    plt.savefig(fig_fld / f"{name}.png")


fld = pathlib.Path("/home/marco/ownCloud/neherlab/code/pangenome-evo/results/ST131")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

# %%
df = pd.read_csv(df_file, index_col=0)
Js = df.index.to_list()
# %%

tree = Phylo.read(tree_file, "newick")
tree.ladderize()
terminal_len = {b.name: b.branch_length for b in tree.get_terminals()}

# %%
two_cat_isolates = {}
for j in Js:
    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    paths = ut.pangraph_to_path_dict(pan)
    path_cat = ut.path_categories(paths)

    if not len(path_cat) == 2:
        continue

    n1, p1, i1 = path_cat[0]
    n2, p2, i2 = path_cat[1]

    if n2 == 1:
        continue

    two_cat_isolates[j] = (p1, p2, i1, i2)


# %%
sdf = df.loc[two_cat_isolates.keys()]

fgfld = fig_fld / "two_cat_isolates"
fgfld.mkdir(exist_ok=True)
for k in two_cat_isolates:
    # move all figures in one place
    # os.
