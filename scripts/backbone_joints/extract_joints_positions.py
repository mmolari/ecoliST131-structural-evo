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
    """Returns a dictionary isolate -> [beg, end, strand] for the junction
    encompassed by the specified backbone edge."""
    ln, rn = edge.left, edge.right
    pos_l, pos_r = pan.blocks[ln.id].alignment.pos, pan.blocks[rn.id].alignment.pos
    pos_dict = defaultdict(lambda: [None, None, None])

    for aln_key in pos_l:
        iso, occ, strand = aln_key
        b, e = pos_l[aln_key]
        assert occ == 1
        same_strand = ln.strand == strand
        start = b if same_strand else e
        idx = 0 if same_strand else 1
        pos_dict[iso][idx] = start
        pos_dict[iso][2] = same_strand

    for aln_key in pos_r:
        iso, occ, strand = aln_key
        b, e = pos_r[aln_key]
        assert occ == 1
        same_strand = rn.strand == strand
        end = e if same_strand else b
        idx = 1 if same_strand else 0
        pos_dict[iso][idx] = end
        assert pos_dict[iso][2] == same_strand

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
        # convention: left node is the one alphabetically smaller
        if e.left.id > e.right.id:
            e = e.invert()
        assert e.left.id < e.right.id

        edges_pos[e.to_str_id()] = find_edges_pos(pan, e)

    # write to file
    with open(args.positions, "w") as f:
        json.dump(edges_pos, f)
