# %%
# script to extract the full core alignment for the pangenome graph.
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
    parser.add_argument("--aln_out", help="output fasta alignment file", type=str)
    parser.add_argument(
        "--reference_iso", help="reference isolate for core-genome order", type=str
    )
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

    # get strain names
    strains = np.sort(pan.strains())

    # core-alignment order
    ref_path = pan.paths[args.reference_iso]
    ref_order = [(b, s) for b, s in zip(ref_path.block_ids, ref_path.block_strands)]
    ref_order = np.array([(b, s) for (b, s) in ref_order if b in core_blocks])

    assert set([b for b, s in ref_order]) == set(core_blocks)

    # extract alignments
    core_alns = defaultdict(str)
    for bid, strand in ref_order:
        # core block
        cbl = pan.blocks[bid]

        alns, occs = cbl.alignment.generate_alignments()

        # alignment dictionary
        aln_dict = {occ[0]: aln for aln, occ in zip(alns, occs)}

        # optionally invert strand
        if not strand:
            aln_dict = {s: aln.reverse_complement() for s, aln in aln_dict.items()}

        for iso in strains:
            core_alns[iso] += "".join(aln_dict[iso])

    # save compressed alignment
    save_aln(core_alns, fname=args.aln_out, strains=strains)

# %%
