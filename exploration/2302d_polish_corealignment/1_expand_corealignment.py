# %%
import numpy as np
import pypangraph as pp
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns

# %%
pg_file = "../../results/ST131/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pg_file)
# %%

# %%
Nstrains = len(pan.paths)
strains = np.sort(pan.strains())

aln_matrices = {}
cons_L = {}

for b in pan.blocks:

    # only core blocks
    if b.is_duplicated() or (b.frequency() < Nstrains):
        continue
    print(b)

    # get restricted alignment and mutated positions
    ordered_occs = sorted(b.alignment.occs, key=lambda x: x[0])
    assert np.all(np.array([o[0] for o in ordered_occs]) == strains)
    pos, snps, occs = b.alignment.extract_nongap_SNPs(which=ordered_occs)
    mask = np.all(np.isin(snps, list("ACGT")), axis=0)
    pos = pos[mask]
    snps = snps[:, mask]

    # save as sparse matrix
    L = len(b.sequence)
    cons_L[b.id] = L
    M = np.zeros((Nstrains, L), dtype=bool)
    for s, p in zip(snps.T, pos):
        M[:, p] = s != b.sequence[p]
    M = sp.coo_matrix(M, dtype=bool)

    aln_matrices[b.id] = M
# %%
core_blocks = list(aln_matrices.keys())
cmap = plt.get_cmap("rainbow")
i, T = 0, len(core_blocks)


def next_color():
    global i
    i += 1
    i %= T
    return cmap(i / (T - 1))


bl_colors = defaultdict(next_color)


fig, axs = plt.subplots(Nstrains, 1, figsize=(6.5, Nstrains * 1.3), sharex=True)

mut_count = []

Ltot = sum(cons_L.values())
step = 1000
bins = np.arange(0, Ltot + step, step)

for n, s in enumerate(strains):
    ax = axs[n]

    p = pan.paths[s]

    blocks = [(b, o) for b, o in zip(p.block_ids, p.block_strands) if b in core_blocks]

    aln = []
    idxs = []
    L = 0
    for b, o in blocks:
        row = sp.csr_matrix(aln_matrices[b]).getrow(n).toarray()[0]
        if not o:
            row = row[::-1]
        aln += list(row)
        idxs += list(np.argwhere(row).flatten() + L)
        L += len(row)

    Ls = np.cumsum([cons_L[b] for b, o in blocks])
    for (b, o), l in zip(blocks, Ls):
        ax.axvline(l, c=bl_colors[b], alpha=0.1, zorder=-1)
    ax.hist(idxs, bins=250, color="k")
    ax.text(0.01, 0.8, s, transform=ax.transAxes)

    ct, bn = np.histogram(idxs, bins=bins)
    mut_count += list(ct)

axs[-1].set_xlabel("core genome projection (bp)")
plt.tight_layout()
plt.savefig("figs/core_genome_density.png", facecolor="w")
plt.show()
# break
# %%
plt.figure(figsize=(10, 4))
plt.hist(mut_count, bins=np.arange(max(mut_count) + 2) - 0.5)
plt.yscale("log")
plt.xlabel("n. SNPs per 1kbp on core genome alignment")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("figs/mut_histogram.png", facecolor="w")
plt.show()
# %%
