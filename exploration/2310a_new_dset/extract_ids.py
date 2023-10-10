# %%
import pathlib
from Bio import SeqIO
import pandas as pd
import numpy as np


# %%
def decide_type(descr, L):
    if "chromosome" in descr.lower():
        return "chromosome"
    elif "plasmid" in descr.lower():
        return "plasmid"
    elif L > 4e6:
        return "chromosome"
    elif L < 1e6:
        return "plasmid"
    else:
        return "unknown"


folder = pathlib.Path("../../config/new_dset/n225-subset-genomes")
dset = []
for fname in folder.glob("*.fna"):
    with open(fname, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            L = len(record.seq)
            descr = record.description
            dset.append(
                {
                    "id": record.id,
                    "asmbl": fname.stem,
                    "len": L,
                    "type": decide_type(descr, L),
                    "description": descr,
                }
            )
dset = pd.DataFrame(dset).set_index(["asmbl", "id"])
dset.to_csv("data/new_dset_accs.csv")

# %%

dset.groupby("asmbl")["len"].count().describe()

# %%

dset["type"].value_counts()

# %%
mask = dset["type"] == "chromosome"
chr_acc = dset[mask].index.get_level_values("id").unique().to_numpy()
chr_acc.tofile("data/new_dset_chr_accs.txt", sep="\n", format="%s")


# %%
