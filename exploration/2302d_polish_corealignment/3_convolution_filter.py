# %%

import pypangraph as pp
import numpy as np
import matplotlib.pyplot as plt
from Bio import Align, SeqIO
import scipy.sparse as sps

# %%

# input:
# pangraph
# guide strain (optional)
# window size
# threshold number of mutations
# (remove all mutations in the threshold window)
# hypothesis: few mutations
# output:
# - diagnostic plot
# - restricted alignment
# - restricted and polished alignment
# %%


def block_to_nogap_aln_matrix(b, occs_order):

    aln, occs = b.alignment.generate_alignments(which=occs_order)

    A = np.vstack([np.array(list(a)) for a in aln])

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


def create_alignment_matrices(pan, block_order):

    # ordered list of strains for assert
    strains = np.sort(pan.strains())

    # create containers
    consensus_matrix = []
    aln_matrix = []
    aln_matrix_idxs = []

    # index of current alignment in core-genome alignment
    left_idx = 0

    for bl in block_order:

        b = pan.blocks[bl]

        print(b)

        # order by strain name
        ordered_occs = sorted(b.alignment.occs, key=lambda x: x[0])
        # check that all strains are there
        assert np.all(np.array([o[0] for o in ordered_occs]) == strains)

        # no-gap alignment matrix
        A = block_to_nogap_aln_matrix(b, ordered_occs)
        lA = A.shape[1]
        idxs = np.arange(lA, dtype=int) + left_idx

        # sparse consensus matrix
        cons = consensus(A)
        S = sps.coo_matrix(A != cons)

        # keep only non-consensus
        is_cons = np.all(A[0, :] == A, axis=0)
        A = A[:, ~is_cons]
        idxs = idxs[~is_cons]

        # append
        aln_matrix.append(A)
        aln_matrix_idxs.append(idxs)
        consensus_matrix.append(S)

        # increase index
        left_idx += lA

    # concatenate
    consensus_matrix = sps.hstack(consensus_matrix)
    aln_matrix = np.hstack(aln_matrix)
    aln_matrix_idxs = np.hstack(aln_matrix_idxs)

    return {
        "consensus_matrix": consensus_matrix,
        "aln_matrix": aln_matrix,
        "aln_matrix_idxs": aln_matrix_idxs,
    }


def aln_matrix_to_fasta(A, strains, filename):
    recs = []
    for n, s in enumerate(strains):
        seq = Align.Seq("".join(list(A[n, :])))
        rec = Align.SeqRecord(seq, id=s, description="")
        recs.append(rec)
    SeqIO.write(recs, filename, "fasta")


def filter_out_idxs(S, window, threshold):
    N, L = S.shape
    remove_idxs = []
    for n in range(N):
        row = S.getrow(n)
        assert np.all(row.data)
        idxs = np.sort(row.indices)

        for i in idxs:
            d1 = (idxs - i) % L
            d2 = (i - idxs) % L
            ds = np.min(np.vstack([d1, d2]), axis=0)

            close_mask = ds <= window
            n_neigh = np.sum(close_mask)
            if n_neigh > threshold:
                remove_idxs += list(idxs[close_mask])
    remove_idxs = np.unique(remove_idxs)
    return remove_idxs


def diagnostic_plot(
    Aidxs,
    remove_idxs,
    L,
    window,
    threshold,
    filename=None,
    show=False,
):

    kept_idxs = list(set(Aidxs) - set(remove_idxs))

    bins = np.arange(L + 10001, step=10000)

    kwargs = {
        "bins": bins,
        "cumulative": True,
        "density": True,
        # "alpha" : 0.5,
        "histtype": "step",
    }

    fig, axs = plt.subplots(3, 1, sharex=False, figsize=(10, 8))

    ax = axs[0]
    window_bins = np.arange(0, L + window + 1, step=window)
    ct1, _ = np.histogram(Aidxs, bins=window_bins)
    ct2, _ = np.histogram(kept_idxs, bins=window_bins)
    ctbins = np.arange(np.max(np.hstack([ct1, ct2])) + 3) - 0.5

    ax.hist(ct1, bins=ctbins, alpha=0.4)
    ax.hist(ct2, bins=ctbins, alpha=0.4)
    ax.set_xlim(left=0)
    ax.axvline(threshold, c="k", ls=":", label="threshold")
    ax.set_yscale("log")
    ax.set_xscale("symlog")
    ax.set_title(
        f"window size = {window} bp, threshold > {threshold} neighbouring SNPs"
    )
    ax.set_xlabel("n. core-genome alignment polymorphic positions per window")
    ax.set_ylabel("n. positions")

    ax = axs[1]
    ax.hist(Aidxs, bins=bins, label="pre-filter")
    ax.hist(remove_idxs, bins=bins, label="removed")
    # ax.hist(kept_idxs, bins=bins, label="post-filter", histtype="step")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlabel("core genome alignment")
    ax.set_ylabel("SNPs per 10kbp")
    ax.set_title(
        f"core alignment size before / after filtering = {len(Aidxs)} / {len(kept_idxs)}"
    )

    ax = axs[2]
    ax.hist(Aidxs, label="pre-filtering", **kwargs)
    ax.hist(remove_idxs, label="removed", **kwargs)
    ax.hist(kept_idxs, label="post-filter", **kwargs)
    ax.legend(loc="upper left")
    ax.set_xlabel("core genome alignment")
    ax.set_ylabel("cumul. distr. of SNPs")

    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename)
    if show:
        plt.show()
    else:
        plt.close(fig)


# %%

pg_file = "../../results/ST131/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pg_file)

Nstrains = len(pan.paths)
strains = np.sort(pan.strains())

# ordered list of blocks
guide_strain = "NZ_CP014522"
coreblock_ids = [
    b.id for b in pan.blocks if (not b.is_duplicated()) and (b.frequency() == Nstrains)
]
block_order = [b for b in pan.to_paths_dict()[guide_strain] if b in coreblock_ids]

# get alignment matrices
aln_Ms = create_alignment_matrices(pan, block_order)

# %%

S = aln_Ms["consensus_matrix"]
A = aln_Ms["aln_matrix"]
Aidxs = aln_Ms["aln_matrix_idxs"]
N, L = S.shape

# %%
window = 1000
threshold = 3

remove_idxs = filter_out_idxs(S, window, threshold)
# %%

diagnostic_plot(
    Aidxs,
    remove_idxs,
    L,
    window,
    threshold,
    filename=None,
    show=True,
)

# %%

keep = ~np.isin(Aidxs, remove_idxs)
A_polished = A[:, keep]
# %%

polished_aln_filename = "aln.fa"
aln_matrix_to_fasta(A_polished, strains, polished_aln_filename)
# %%
