# %%
from Bio import SeqIO
import pandas as pd

fname = "data/hochhauser.csv"
df = pd.read_csv(fname, index_col=0)
df

# %%
upstr = df.iloc[:, 0].to_dict()
downstr = df.iloc[:, 3].to_dict()
locs = list(set(upstr.values()) | set(downstr.values()))

# %%
fname = "data/IMG_2687453259/IMG_Data/2687453259/2687453259.genes.fna"
records = []
with open(fname) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if int(record.id) in locs:
            records.append(record)

# %%
svname = "data/anchor_genes.fa"
with open(svname, "w") as handle:
    SeqIO.write(records, handle, "fasta")

# %%
