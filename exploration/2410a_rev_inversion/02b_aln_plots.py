# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib
from collections import Counter
import pypangraph as pp
from Bio import AlignIO, Phylo
import gzip

# look into inversion that happened multiple times
# 115 kbps
# isolates:
# NZ_CP107114.1
# NZ_CP107184.1
# NZ_CP107182.1

fld = pathlib.Path("../../results/ST131_ABC")
res_fld = pathlib.Path("res")
res_fld.mkdir(exist_ok=True)

fig_fld = pathlib.Path("figs/f02")
fig_fld.mkdir(exist_ok=True, parents=True)

colors_file = fld / "pangraph/coresynt-asm20-100-5/blocks.csv"
mergers_file = fld / "pangraph/coresynt-asm20-100-5/mergers.csv"

cl_df = pd.read_csv(colors_file, index_col=0)
mg_df = pd.read_csv(mergers_file, index_col=0)["0"]

pangraph_file = fld / "pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
isolates = pan.strains()

tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, format="newick")

iso_A = "NZ_CP107114.1"
iso_B = "NZ_CP107184.1"
iso_C = "NZ_CP107182.1"

# %%

pattern = "full"
aln_file = res_fld / f"{pattern}_alignment_ungapped.fa.gz"
with gzip.open(aln_file, "rt") as f:
    aln = AlignIO.read(f, format="fasta")

pos_file = res_fld / f"{pattern}_alignment_block_pos.csv"
Bpos = pd.read_csv(pos_file)

starts = Bpos.groupby("merger")["start"].min()
ends = Bpos.groupby("merger")["end"].max()
SE = pd.DataFrame({"start": starts, "end": ends}).sort_values("start")
# %%
M = np.array(aln)


# %%
# alignment consensus
def consensus(aln_col):
    c = Counter(aln_col)
    return c.most_common()[0][0]


cons = np.apply_along_axis(consensus, 0, M)
# %%
B = M != cons
N = M.shape[0]
# sum chunks of 1kbp

chunk_size = 1000
X = B.shape[1] // chunk_size
Bsum = np.zeros((N, X))
for x in range(X):
    s, e = x * chunk_size, (x + 1) * chunk_size
    S = np.sum(B[:, s:e], axis=1)
    Bsum[:, x] = S
# %%

strains = [x.name for x in tree.get_terminals()]
aln_isos = [x.id for x in aln]

# %%
fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharey=True, sharex=True)
for n, iso in enumerate([iso_A, iso_B, iso_C]):
    ax = axs[n]
    idx = aln_isos.index(iso)

    x = np.arange(X) * chunk_size
    y = Bsum[idx] / chunk_size
    ax.plot(x, y)
    # ax.plot(x, np.zeros_like(x), color="w", lw=2)

    for i, row in SE.iterrows():
        s = row["start"]
        if s == 0:
            continue
        ax.axvline(s, color="k")

    ax.set_xlabel(f"core alignment {iso} (bp)")
    ax.set_ylabel("nonconsensus density")
    ax.set_ylim(bottom=0)
    ax.set_xlim(0, B.shape[1])
plt.tight_layout()
plt.savefig(fig_fld / "aln_nonconsensus.png", dpi=300)
plt.show()

# %%
