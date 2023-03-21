# %%

import pypangraph as pp
import numpy as np
import matplotlib.pyplot as plt
from Bio import Align, SeqIO
import scipy.sparse as sps

# %%

pg_file = "../../results/ST131/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pg_file)
# %%


def block_to_nogap_aln_matrix(b, occs_order):

    aln, occs = b.alignment.generate_alignments(which=occs_order)

    recs = []
    for a, o in zip(aln, occs):
        rec = Align.SeqRecord(a, name=o[0])
        recs.append(rec)
    A = Align.MultipleSeqAlignment(recs)
    A = np.array(A)

    gap_mask = np.any(A == "-", axis=0)
    nuc_mask = np.all(np.isin(A, list("ACGT")), axis=0)
    mask = nuc_mask & (~gap_mask)
    return A[:, mask]


def consensus(A):
    def line(l):
        lett, ct = np.unique(l, return_counts=True, axis=0)
        return lett[np.argmax(ct)]

    cons = np.apply_along_axis(line, 0, A)
    return cons


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

    # order by strain name
    ordered_occs = sorted(b.alignment.occs, key=lambda x: x[0])
    assert np.all(np.array([o[0] for o in ordered_occs]) == strains)

    A = block_to_nogap_aln_matrix(b, ordered_occs)
    cons = consensus(A)

    cons_L[b.id] = len(cons)

    S = sps.coo_matrix(A != cons)
    aln_matrices[b.id] = S

# %%
guide_strain = "NZ_CP014522"
block_order = [b for b in pan.to_paths_dict()[guide_strain] if b in aln_matrices]
K = sps.hstack([aln_matrices[b] for b in block_order])
L = K.shape[1]
step = 1000
bins = np.arange(0, L + step + 2, step=step)
# %%
mask_binct = np.ones(len(bins) - 1, dtype=bool)
cts = []
for n in range(Nstrains):
    row_idx = K.getrow(n).indices
    ct, bn, pl = plt.hist(
        row_idx,
        bins=bins,
        histtype="step",
        # density=True,
        # cumulative=True
        alpha=0.5,
    )
    cts += list(ct)
    mask_binct[ct >= 3] = False
plt.yscale("log")
plt.xlabel("core genome alignment position")
plt.ylabel("n. SNPs per kbp")
plt.show()
# %%

plt.hist(cts, bins=np.arange(np.max(cts) + 2) - 0.5)
plt.xscale("symlog")
plt.yscale("log")
plt.xlabel("SNPs per kbp")
plt.ylabel("n. of 1kbp bins")
plt.show()
# %%
plt.hist(K.sum(axis=1), bins=np.logspace(1, 4, 100), cumulative=True)
plt.xscale("log")
plt.xlabel("n. mutations per genome")
plt.ylabel("n. genomes")
plt.show()
# %%
Dthreshold = 1000

dL = []
for n in range(Nstrains):
    r = K.getrow(n)
    deltas = np.diff(np.sort(r.indices))
    dL += list(deltas)

plt.hist(dL, bins=np.logspace(0, 7, 100))
plt.xscale("log")
plt.axvline(Dthreshold, c="k", ls=":")
plt.xlabel("distance between subsequent mutations")
plt.show()
# %%

mask = np.ones(L, dtype=bool)

for n in range(Nstrains):
    idxs = np.sort(K.getrow(n).indices)
    deltas = np.diff(list(idxs) + [idxs[0] + L])
    remove = (deltas < Dthreshold) & (np.roll(deltas, 1) < Dthreshold)
    rem_idxs = np.argwhere(remove).flatten()
    print(deltas[rem_idxs])
    print(deltas[rem_idxs - 1])

    for ri in rem_idxs:
        si = (ri - 1) % len(idxs)
        ei = (ri + 1) % len(idxs)

        if si > ei:
            mask[idxs[si] :] = False
            mask[: idxs[ei]] = False
        else:
            mask[idxs[si] : idxs[ei]] = False

# %%
plt.figure(figsize=(10, 2))
plt.plot(mask, label="mutation distance")
plt.plot(bins[:-1], -mask_binct.astype(int), label="mutations per bin")
plt.yticks([-1, 0, 1], labels=["keep", "exclude", "keep"])
plt.legend()
plt.xlabel("core genome alignment")
plt.tight_layout()
plt.show()
# %%

As = []
SNPs_idxs = []

left_idx = 0

for bl in block_order:

    b = pan.blocks[bl]

    print(b)

    # order by strain name
    ordered_occs = sorted(b.alignment.occs, key=lambda x: x[0])
    assert np.all(np.array([o[0] for o in ordered_occs]) == strains)

    A = block_to_nogap_aln_matrix(b, ordered_occs)

    keep_idxs = np.arange(A.shape[1], dtype=int) + left_idx

    l = A.shape[1]
    sub_mask = mask[left_idx : left_idx + l]
    A = A[:, sub_mask]
    keep_idxs = keep_idxs[sub_mask]

    is_cons = np.all(A[0, :] == A, axis=0)
    A = A[:, ~is_cons]
    keep_idxs = keep_idxs[~is_cons]

    As.append(A)

    SNPs_idxs += list(keep_idxs)
    left_idx += l

As = np.hstack(As)

# %%
recs = []
for n, s in enumerate(strains):
    seq = Align.Seq("".join(list(As[n, :])))
    rec = Align.SeqRecord(seq, id=s, description="")
    recs.append(rec)


# %%
SeqIO.write(recs, "reduced_aln.fa", "fasta")
# %%
all_snps_idxs = np.unique(K.nonzero()[1])
plt.hist(SNPs_idxs, bins=1000, cumulative=True, density=True, label="after refining")
plt.hist(all_snps_idxs, bins=1000, cumulative=True, density=True,
    # histtype="step",
    label="before refining",
    alpha=0.5,
)
plt.xlabel("core genome alignment")
plt.ylabel("cumulative mutation distribution")
plt.legend()
# plt.plot(mask)
plt.show()
# %%
