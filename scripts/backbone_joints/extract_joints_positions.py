import argparse
import json
import pandas as pd

import pypangraph as pp
import utils as ut

from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Given a pangraph and a list of backbone edges, extracts the
        position of backbone joints beginning and end.
        """
    )
    parser.add_argument("--pangraph", type=str, required=True)
    parser.add_argument("--edges", type=str, required=True)
    parser.add_argument("--positions", type=str, required=True)
    return parser.parse_args()


def find_edges_pos(pan, edge):
    """Returns a dictionary isolate -> [left beg, left end, right beg, right end, strand]
    for the junction encompassed by the specified backbone edge.
    Note that lb < le and rb < re unless one wraps around the genome.
    The area of interest always sits between lb -> re moving in increasing order,
    even if it wraps around the genome.
    """
    ln, rn = edge.left, edge.right
    pos_l, pos_r = pan.blocks[ln.id].alignment.pos, pan.blocks[rn.id].alignment.pos
    pos_dict = defaultdict(lambda: [None, None, None, None, None])

    for aln_key in pos_l:
        iso, occ, strand = aln_key
        beg, end = pos_l[aln_key]
        assert occ == 1
        same_strand = ln.strand == strand
        bidx, eidx = (0, 1) if same_strand else (2, 3)
        pos_dict[iso][bidx] = beg
        pos_dict[iso][eidx] = end
        pos_dict[iso][4] = same_strand

    for aln_key in pos_r:
        iso, occ, strand = aln_key
        beg, end = pos_r[aln_key]
        assert occ == 1
        same_strand = rn.strand == strand
        bidx, eidx = (2, 3) if same_strand else (0, 1)
        pos_dict[iso][bidx] = beg
        pos_dict[iso][eidx] = end
        assert pos_dict[iso][4] == same_strand

    return dict(pos_dict)


if __name__ == "__main__":
    args = parse_args()

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)
    bdf = pan.to_blockstats_df()
    is_core = bdf["core"].to_dict()
    block_len = bdf["len"].to_dict()
    strains = pan.strains()
    N_iso = len(strains)

    # backbone blocks and edges
    df = pd.read_csv(args.edges)
    df = df[df["count"] == N_iso]
    backbone_edges = [ut.Edge.from_str_id(e) for e in df["edge"]]

    # find edges positions
    edges_pos = {}
    for e in backbone_edges:
        edges_pos[e.to_str_id()] = find_edges_pos(pan, e)

    # write to file
    with open(args.positions, "w") as f:
        json.dump(edges_pos, f)
