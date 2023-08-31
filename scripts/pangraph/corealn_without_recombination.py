# %%
# script to extract the reduced core genome alignment (only SNPs),
# but with a filter on the maximum number of mutations per kbp, to
# remove recombination spots.

import argparse
import pypangraph as pp
import numpy as np
import json
import scipy.sparse as sps
from Bio import SeqIO, Seq


# %%
def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Produces a polished core genome alignment where regions with high SNPs density are
    removed. This is done to filter out recombination spots.
    """
    )
    parser.add_argument("--pangraph", help="json pangraph file", type=str)
    parser.add_argument("--guide_strain", help="determines core genome order", type=str)
    parser.add_argument("--window", help="window size (bp)", type=int)
    parser.add_argument("--max_nsnps", help="max n. of SNPs in the window", type=int)
    parser.add_argument("--fasta_aln", help="output fasta alignment file", type=str)
    parser.add_argument("--info_size", help="output json information file", type=str)
    parser.add_argument("--info_idxs", help="save idxs of reduced alignment", type=str)
    args = parser.parse_args()
    return args


def block_id_order(pan, guide_strain):
    """Given a guide strain and a pangraph returns core blocks in the order
    they appear in the strain."""
    N = len(pan.paths)
    is_core = lambda b: (not b.is_duplicated()) and (b.frequency() == N)
    coreblock_ids = [b.id for b in pan.blocks if is_core(b)]
    block_order = [b for b in pan.to_paths_dict()[guide_strain] if b in coreblock_ids]
    return block_order


def block_to_aln(b, occs_order):
    """given a block and an order of occurrences returns a reduced alignment
    matrix where gaps have been stripped away."""
    aln, occs = b.alignment.generate_alignments(which=occs_order)
    A = np.vstack([np.array(list(a)) for a in aln])
    gap_mask = np.any(A == "-", axis=0)
    nuc_mask = np.all(np.isin(A, list("ACGT")), axis=0)
    mask = nuc_mask & (~gap_mask)
    return A[:, mask]


def consensus(A):
    """given an alignment matrix returns the consensus"""

    def site_consensus(l):
        lett, ct = np.unique(l, return_counts=True, axis=0)
        return lett[np.argmax(ct)]

    return np.apply_along_axis(site_consensus, 0, A)


def create_alignment_matrices(pan, block_order):
    # alphabetically ordered list of strains
    strains = np.sort(pan.strains())

    # create containers
    consensus_matrix = []
    aln_matrix = []
    aln_matrix_idxs = []

    # index of current alignment in core-genome alignment
    left_idx = 0

    for bl in block_order:
        b = pan.blocks[bl]

        print("proccessing", b)

        # order by strain name
        ordered_occs = sorted(b.alignment.occs, key=lambda x: x[0])
        # check that all strains are there
        assert np.all(np.array([o[0] for o in ordered_occs]) == strains)

        # no-gap alignment matrix
        A = block_to_aln(b, ordered_occs)
        lA = A.shape[1]
        idxs = np.arange(lA, dtype=int) + left_idx

        # sparse consensus matrix
        cons = consensus(A)
        S = sps.coo_matrix(A != cons)

        # keep only non-consensus columns
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
    # reduced alignment
    aln_matrix = np.hstack(aln_matrix)
    # indices of SNPs in the full core-genome alignment
    aln_matrix_idxs = np.hstack(aln_matrix_idxs)
    # boolean sparse matrix of SNP (1) and consensus (0)
    consensus_matrix = sps.hstack(consensus_matrix)

    return {
        "consensus_matrix": consensus_matrix,
        "aln_matrix": aln_matrix,
        "aln_matrix_idxs": aln_matrix_idxs,
        "strain_order": strains,
    }


def filter_out_idxs(S, window, max_nsnps):
    """Filters out sites where mutations are too clustered on a
    specific strain. This is defined using a window size and an upper
    threshold n. of SNPs per window"""
    N, L = S.shape
    remove_idxs = []
    for n in range(N):
        row = S.getrow(n)
        assert np.all(row.data)
        idxs = np.sort(row.indices)

        for i in idxs:
            # ds -> vector of distances of SNP positions to window center
            # with periodic boundary conditions
            d1 = (idxs - i) % L
            d2 = (i - idxs) % L
            ds = np.min(np.vstack([d1, d2]), axis=0)

            # boolean array of which SNPs are within the window
            close_mask = ds <= window
            n_neigh = np.sum(close_mask)
            if n_neigh > max_nsnps:
                remove_idxs += list(idxs[close_mask])
    # return list of indices to remove
    remove_idxs = np.unique(remove_idxs)
    return remove_idxs


def save_aln(A, strains, fname):
    """Saves an alignment matrix in a fasta format"""
    recs = []
    for n, s in enumerate(strains):
        seq = Seq.Seq("".join(list(A[n, :])))
        rec = SeqIO.SeqRecord(seq, id=s, description="")
        recs.append(rec)
    SeqIO.write(recs, fname, format="fasta")


def save_summary_info(S, remove_idxs, window, nsnps_max, fname_size, fname_idxs):
    """Saves infos on the size of the core genome alignment before and
    after filtering."""

    N, L = S.shape
    Aidxs = np.unique(S.col)
    keep_idxs = list(set(Aidxs) - set(remove_idxs))

    # n. of total remaining sites
    mask = np.ones(L, dtype=bool)
    for i in remove_idxs:
        s, e = (i - window) % L, (i + window + 1) % L
        if s < e:
            mask[s:e] = False
        else:
            mask[s:] = False
            mask[:e] = False
    L_filtered = np.sum(mask)

    info = {
        "core aln size": L,
        "core aln snps": len(Aidxs),
        "core aln consensus": L - len(Aidxs),
        "polished aln size": L_filtered,
        "polished aln snps": len(keep_idxs),
        "polished aln consensus": L_filtered - len(keep_idxs),
        "window size": window,
        "nsnps max": nsnps_max,
    }

    # make json serializable
    info = {k: int(v) for k, v in info.items()}

    with open(fname_size, "w") as f:
        json.dump(info, f, indent=2)

    # save idxs
    info = {
        "idxs_all": Aidxs.tolist(),
        "idxs_removed": remove_idxs.tolist(),
    }
    with open(fname_idxs, "w") as f:
        json.dump(info, f, indent=2)


if __name__ == "__main__":
    # parse arguments
    args = parse_args()

    # load pangenome graph
    pan = pp.Pangraph.load_json(args.pangraph)

    # ordered list of blocks
    block_order = block_id_order(pan, args.guide_strain)

    # get alignment matrices
    aln_Ms = create_alignment_matrices(pan, block_order)

    S = aln_Ms["consensus_matrix"]
    A = aln_Ms["aln_matrix"]
    Aidxs = aln_Ms["aln_matrix_idxs"]
    strains = aln_Ms["strain_order"]
    N, L = S.shape

    # filter out highly mutated (prbably recombined) spots
    w, thr = args.window, args.max_nsnps
    remove_idxs = filter_out_idxs(S, window=w, max_nsnps=thr)

    # remove recombination islands
    keep = ~np.isin(Aidxs, remove_idxs)
    A_polished = A[:, keep]

    # save polished alignment
    save_aln(A_polished, strains=strains, fname=args.fasta_aln)

    # save summary info
    save_summary_info(
        S,
        remove_idxs,
        window=w,
        nsnps_max=thr,
        fname_size=args.info_size,
        fname_idxs=args.info_idxs,
    )
