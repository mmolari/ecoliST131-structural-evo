# %%
import pathlib
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import Phylo

tree_file = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, format="newick")
tree.ladderize()
str_ord = [x.name for x in tree.get_terminals()]

df_file = "../../results/ST131_ABC/distances/summary-asm20-100-5.csv"
adf = pd.read_csv(df_file)

# %%

Phylo.draw(tree.root[1][1][1][1], label_func=lambda x: "")

A = [t.name for t in tree.root[0].get_terminals()]
B = [t.name for t in tree.root[1][0].get_terminals()]
C = [t.name for t in tree.root[1][1].get_terminals()]
# %%
for clade in [A, B, C]:
    print(len(clade))
    mask = adf["si"].isin(clade) & adf["sj"].isin(clade) & (adf["si"] > adf["sj"])
    values = adf[mask]["core_div_filtered"]
    print(values.mean() * 1e6)
mask = adf["si"] > adf["sj"]
values = adf[mask]["core_div_filtered"]
print(values.mean() * 1e6)

# %%
