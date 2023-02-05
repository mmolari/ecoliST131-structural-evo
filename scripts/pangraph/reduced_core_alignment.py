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

    return {"gaps": n_gaps, "consensus": n_consensus, "SNPs_matrix": M}


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
    }

    # get strain names
    strains = np.sort(pan.strains())

    # extract alignments
    n_gaps = 0
    n_consensus = 0
    n_snps = 0
    core_alns = defaultdict(str)
    for cbl_name in core_blocks:
        # core block
        cbl = pan.blocks[cbl_name]

        alns, occs = cbl.alignment.generate_alignments()

        # alignment dictionary
        aln_dict = {occ[0]: aln for aln, occ in zip(alns, occs)}

        # reduce alignment
        red_aln = reduce_aln(aln_dict, strains)
        n_gaps += red_aln["gaps"]
        n_consensus += red_aln["consensus"]
        M = red_aln["SNPs_matrix"]
        n_snps += M.shape[1]
        for n, s in enumerate(strains):
            core_alns[s] += "".join(M[n])

    # add and save info
    info["n. gaps"] = int(n_gaps)
    info["n. consensus"] = int(n_consensus)
    info["n. snps"] = int(n_snps)
    print(info)
    with open(args.info, "w") as f:
        json.dump(info, f)

    # save compressed alignment
    save_aln(core_alns, fname=args.fasta_aln, strains=strains)

# %%
