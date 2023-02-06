# %%
# script to extract the reduced core alignment for the pangenome graph.
# this can be used to build a phylogenetic tree / distance matrix

import argparse
import pypangraph as pp
import numpy as np
import json
from Bio import SeqIO, Seq
from collections import defaultdict

# %%
def parse_args():
    parser = argparse.ArgumentParser(description="extract core blocks alignments")
    parser.add_argument("--pangraph", help="json pangraph file", type=str)
    parser.add_argument("--fasta_aln", help="output fasta alignment file", type=str)
    parser.add_argument("--info", help="output json information file", type=str)
    args = parser.parse_args()
    return args


def save_aln(aln_dict, fname, strains):
    """Given a dictionary {strain -> sequence} saves the alignment in a fasta file
    with order specified by the `strains` argument."""
    recs = []
    for strain in strains:
        seq = Seq.Seq(aln_dict[strain])
        rec = SeqIO.SeqRecord(seq, id=strain, name="", description="")
        recs.append(rec)
    SeqIO.write(recs, fname, format="fasta")


def reduce_aln(aln_dict, strains):
    """Given an alignment dictionary {strain -> alignment sequence}
    it restrict the alignment to all columns with relevant phylogenetic
    information, excluding gaps and special characters."""

    L = len(aln_dict[strains[0]])
    N = len(strains)
    M = np.empty((N, L), dtype=str)

    for n, s in enumerate(strains):
        M[n] = list(aln_dict[s])

    # remove gaps
    has_gap = np.any(M == "-", axis=0)
    n_gaps = has_gap.sum()
    M = M[:, ~has_gap]

    # remove consensus
    is_consensus = np.all(M == M[0], axis=0)
    n_consensus = is_consensus.sum()
    M = M[:, ~is_consensus]

    # remove extra characters
    is_extra = ~np.all(np.isin(M, list("acgtACGT")), axis=0)
    n_extra = is_extra.sum()
    M = M[:, ~is_extra]

    return {
        "gaps": int(n_gaps),
        "consensus": int(n_consensus),
        "n. special characters": int(n_extra),
        "n. snps": int(M.shape[1]),
        "SNPs_matrix": M,
    }


# %%
if __name__ == "__main__":

    # parse arguments
    args = parse_args()

    # load pangenome graph
    pan = pp.Pangraph.load_json(args.pangraph)

    # extract list of core blocks
    df = pan.to_blockstats_df()
    dfc = df[df["core"]].sort_values("len", ascending=False)
    core_blocks = dfc.index.to_numpy()

    # info dictionary
    info = {
        "n. core blocks": len(dfc),
        "cumulative core length": int(dfc["len"].sum()),
        "n. gaps": 0,
        "n. consensus": 0,
        "n. snps": 0,
        "n. special characters": 0,
    }

    # get strain names
    strains = np.sort(pan.strains())

    # extract alignments
    core_alns = defaultdict(str)
    for cbl_name in core_blocks:
        # core block
        cbl = pan.blocks[cbl_name]

        alns, occs = cbl.alignment.generate_alignments()

        # alignment dictionary
        aln_dict = {occ[0]: aln for aln, occ in zip(alns, occs)}

        # reduce alignment
        red_aln = reduce_aln(aln_dict, strains)

        # append to compressed core alignment
        M = red_aln["SNPs_matrix"]
        for n, s in enumerate(strains):
            core_alns[s] += "".join(M[n])

        # update information
        info["n. gaps"] += red_aln["gaps"]
        info["n. consensus"] += red_aln["consensus"]
        info["n. snps"] += red_aln["n. snps"]
        info["n. special characters"] += red_aln["n. special characters"]

    # add and save info
    print(info)
    with open(args.info, "w") as f:
        json.dump(info, f, indent=2)

    # save compressed alignment
    save_aln(core_alns, fname=args.fasta_aln, strains=strains)

# %%
