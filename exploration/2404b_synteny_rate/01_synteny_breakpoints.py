# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import pathlib

dest = pathlib.Path("data/synt_map")
dest.mkdir(exist_ok=True, parents=True)

ref = "NZ_CP096110.1"
# ref = "NZ_CP124495.1"
qry = "NZ_CP107114.1"
# qry = "NZ_CP107184.1"
# qry = "NZ_CP107182.1"
# qry = "NZ_CP059281.1"
# qry = "NZ_CP059279.1"
# qry = "NZ_OX637961.1"

ref_fname = f"../../data/fa/{ref}.fa"
qry_fname = f"../../data/fa/{qry}.fa"
out_fname = dest / f"{ref}_{qry}.paf"

# map with minimap2
cmd = f"minimap2 -x asm5 {ref_fname} {qry_fname} > {out_fname}"
subprocess.run(cmd, shell=True)

# %%
# load paf results
cols = [
    "qry",
    "qry_len",
    "qry_st",
    "qry_en",
    "strand",
    "ref",
    "ref_len",
    "ref_st",
    "ref_en",
    "n_match",
    "aln_len",
    "mapq",
]
df = pd.read_csv(out_fname, sep="\t", header=None, usecols=range(12), names=cols)
df

# %%

# dotplot
fig, ax = plt.subplots(figsize=(10, 10))
for idx, row in df.iterrows():
    if row["mapq"] < 50:
        continue
    x = [row["qry_st"], row["qry_en"]]
    y = [row["ref_st"], row["ref_en"]]
    if row["strand"] == "-":
        x = x[::-1]
    ax.plot(x, y, lw=1, c="k", marker=".", markersize=2)
    ax.axvline(row["qry_st"], lw=0.5, c="C0", ls="--", alpha=0.5)
    ax.axhline(row["ref_st"], lw=0.5, c="C0", ls="--", alpha=0.5)

ax.set_xlabel(f"{qry} position")
ax.set_ylabel(f"{ref} position")
ax.set_xlim(0, df["qry_len"].max())
ax.set_ylim(0, df["ref_len"].max())
# plt.grid()
plt.tight_layout()
plt.savefig(dest / f"{ref}_{qry}.png", dpi=150)
plt.show()

# %%
df[df["mapq"] > 50]

# qry = "NZ_CP107114.1"
# qry 2639187	2776104
# ref 4484606	4621524

# qry = "NZ_CP107184.1"
# qry 4438103	4587942
# ref 4471684	4621524

# qry = "NZ_CP107182.1"
# qry 328587	478433
# ref 4471677	4621524

# %%
