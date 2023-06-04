# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pathlib
from Bio import Phylo, SeqIO
import pypangraph as pp

import utils as ut

from collections import defaultdict
from itertools import combinations

# %%

svfld = pathlib.Path("figs")
svfld.mkdir(exist_ok=True)


def svfig(name):
    plt.savefig(svfld / name, dpi=300, facecolor="white", bbox_inches="tight")


pair_1 = ["NZ_CP084679", "NZ_CP084678"]
pair_2 = ["NZ_CP035516", "NZ_CP035477"]
pairs = pair_1 + pair_2


# %%
prefix = "../../results/ST131/pangraph"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)
# %%

df = pan.to_blockcount_df()

df = df.loc[pairs]

# remove columns which are all zeros
df = df.loc[:, (df != 0).any(axis=0)]

# select columns that are all ones
anchor_genes = df.loc[:, (df == 1).all(axis=0)].columns.to_numpy()

# select duplicated genes: columns that have values >1 at least once
dupl_genes = df.loc[:, (df > 1).any(axis=0)].columns.to_numpy()

# select accessory genes: columns that are neither anchor nor duplicated
acc_genes = df.loc[
    :, ~df.columns.isin(np.concatenate([anchor_genes, dupl_genes]))
].columns.to_numpy()

# %%

Ls = pan.to_blockstats_df()["len"].to_dict()

# %%

max_anchor = anchor_genes[np.argmax([Ls[a] for a in anchor_genes])]

# %%

cmap = plt.get_cmap("jet")
cdict = defaultdict(lambda: cmap(np.random.rand()))

fig, ax = plt.subplots(1, 1, figsize=(50, 3))

for n, strain in enumerate(pairs):
    x = 0
    p = pan.paths[strain]
    bls = list(p.block_ids)
    i = bls.index(max_anchor)
    if not p.block_strands[i]:
        bls = bls[::-1]
    i = bls.index(max_anchor)
    bls = bls[i:] + bls[:i]

    for b in bls:
        if b in anchor_genes:
            continue
        if b in dupl_genes:
            continue
        l = Ls[b]
        c = cdict[b]
        plt.plot([x, x + l], [n, n], color=c)
        x += l

plt.tight_layout()
svfig("test.png")
plt.show()

# %%
mpair = pair_2
cmd = f"""
pangraph marginalize -s {mpair[0]},{mpair[1]} {pangraph_file} \
    > marginal_2.json
pangraph export -nd -p export_2 -ell 0 marginal_2.json
"""
import os

os.system(cmd)

# %%
# pair 1
# weird_pairs = ["GDCLUOAGIC", "BFCMHGYFHQ"]
# weird_pairs = ["NUCBMDMWWW", "LRNPTTDYZW"]
# weird_pairs = ["DLRPOLYIJF", "XRNULSWROH"]
# weird_pairs = ["LTDLVBKPWI", "OOAGQUOHVB"]

# pair 2
# ~ 300 bps with ~15 % divergence and indels
# weird_pairs = ["HWQDTXMAGV", "CHAMZVRLJX"]
weird_pairs = ["HWQDTXMAGV", "CHAMZVRLJX"]

fa = list(SeqIO.parse("export/export_2.fa", format="fasta"))
fa = {s.name: s for s in fa}

# %%
SeqIO.write(
    [fa[s] for s in weird_pairs],
    f"aln_test/{'_'.join(weird_pairs)}.fa",
    format="fasta",
)

# %%
