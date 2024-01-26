# %%
import pandas as pd
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import pathlib

# %%

fld = pathlib.Path("../../results/ST131_ABC/")

data_fld = pathlib.Path("data")
data_fld.mkdir(parents=True, exist_ok=True)

# fig_fld = pathlib.Path("figs/f00")
# fig_fld.mkdir(parents=True, exist_ok=True)

tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]

integron_summary_file = fld / "annotations/integron_finder/integron_annotations.tsv"
df = pd.read_csv(integron_summary_file, sep="\t")
df["integron_n"] = df["ID_integron"].str.replace("integron_", "").astype(int)
df

# %%
mask = df["type"] == "CALIN"
df[mask]["annotation"].value_counts()
# %%
mask = df["type"] == "complete"
df[mask]["annotation"].value_counts()

# %%


beg = df.groupby(["ID_replicon", "integron_n"])["pos_beg"].min()
end = df.groupby(["ID_replicon", "integron_n"])["pos_end"].max()
tp = df.groupby(["ID_replicon", "integron_n"])["type"].first()

idf = pd.DataFrame({"beg": beg, "end": end, "type": tp})

for idx, row in idf.iterrows():
    if np.abs(row["beg"] - row["end"]) > 1e6:
        raise ValueError(f"Integron {idx} is too big")

idf["len"] = idf["end"] - idf["beg"] + 1
idf = idf.reset_index()
idf.rename(columns={"ID_replicon": "iso"}, inplace=True)
idf["idx"] = idf["iso"] + "|" + idf["integron_n"].astype(str) + "|" + idf["type"]
idf.set_index("idx", inplace=True)
idf.to_csv(data_fld / "integron_regions.csv")
idf
# %%
