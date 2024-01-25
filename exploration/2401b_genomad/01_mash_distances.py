# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import subprocess
from Bio import SeqIO

data_fld = pathlib.Path("data/f01")
data_fld.mkdir(exist_ok=True, parents=True)

fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(exist_ok=True, parents=True)

# %%

# input files
iso1, n1 = "NZ_CP103755.1", 9
# iso1, n1 = "NZ_CP063515.1", 7

iso2 = "NZ_CP049077.2"
# iso2 = "NZ_CP107178.1"
# iso2 = "NZ_CP124455.1"
# iso2 = "NZ_CP063515.1"
f1 = f"../../data/genomad/{iso1}/{iso1}_summary/{iso1}_virus.fna"
f2 = f"../../data/genomad/{iso2}/{iso2}_summary/{iso2}_virus.fna"

# save in single fasta file
seqs_file = data_fld / "prophages.fasta"
cmd = f"cat {f1} {f2} > {seqs_file}"
subprocess.run(cmd, shell=True)

# check sequence length
seq_ls = {}
for seq_record in SeqIO.parse(seqs_file, "fasta"):
    seq_ls[seq_record.id] = len(seq_record.seq)
    print(seq_record.id, len(seq_record.seq))

# run mash
mash_out = data_fld / "mash_tri.tsv"
cmd = f"mash triangle -r {seqs_file} > {mash_out}"
subprocess.run(cmd, shell=True)


# %%

thr = 0.1

# parse mash output
with open(mash_out, "r") as f:
    lines = f.readlines()
N = int(lines[0].strip())
names = []
dist = np.zeros((N, N))
for i in range(N):
    l = lines[i + 1].strip().split("\t")
    name = l[0].replace("provirus_", "")
    iso, pos = name.split("|")
    start, end = pos.split("_")
    name = f"{iso}|{start:0>7}-{end:0>7}"

    names.append(name)
    dist[i, :i] = np.array(l[1:], dtype=float)
    # dist[:i, i] = np.array(l[1:], dtype=float)
    dist[:i, i] = np.nan
    dist[i, i] = np.nan
df = pd.DataFrame(dist, index=names, columns=names)

fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(
    df[df < thr],
    cmap="viridis_r",
    ax=ax,
    cbar_kws={"label": "Mash distance"},
    vmin=0,
    vmax=thr,
)
plt.plot([0, N], [0, N], color="gray")
for i in [0, n1, N]:
    plt.plot([0, i], [i, i], color="gray")
    plt.plot([i, i], [i, N], color="gray")
plt.tight_layout()
plt.savefig(fig_fld / f"mash_tri_{iso1}_{iso2}.png", dpi=300)
plt.show()
# %%
