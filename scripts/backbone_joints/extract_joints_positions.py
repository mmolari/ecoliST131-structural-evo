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
    parser.add_argument("--edge_len_df", type=str, required=True)
    parser.add_argument("--positions", type=str, required=True)
    return parser.parse_args()


def find_edges_pos(pan, edge, isolates):
    """Returns a dictionary isolate -> [left beg, left end, right beg, right end, strand]
    for the junction encompassed by the specified backbone edge.
    Note that lb < le and rb < re unless one wraps around the genome.
    The area of interest always sits between lb -> re moving in increasing order,
    even if it wraps around the genome.
    """
    ln, rn = edge.left, edge.right
    pos_l, pos_r = pan.blocks[ln.id].alignment.pos, pan.blocks[rn.id].alignment.pos
    pos_dict = defaultdict(dict)

    for aln_key in pos_l:
        iso, occ, strand = aln_key
        if iso not in isolates:
            continue
        beg, end = pos_l[aln_key]
        assert occ == 1
        same_strand = ln.strand == strand
        side = "left" if same_strand else "right"
        pos_dict[iso][side + "_beg"] = beg
        pos_dict[iso][side + "_end"] = end
        pos_dict[iso]["strand"] = same_strand

    for aln_key in pos_r:
        iso, occ, strand = aln_key
        if iso not in isolates:
            continue
        beg, end = pos_r[aln_key]
        assert occ == 1
        same_strand = rn.strand == strand
        side = "right" if same_strand else "left"
        pos_dict[iso][side + "_beg"] = beg
        pos_dict[iso][side + "_end"] = end
        assert pos_dict[iso]["strand"] == same_strand

    return {
        iso: (d["left_beg"], d["left_end"], d["right_beg"], d["right_end"], d["strand"])
        for iso, d in pos_dict.items()
    }


if __name__ == "__main__":
    args = parse_args()

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)
    bdf = pan.to_blockstats_df()
    is_core = bdf["core"].to_dict()
    block_len = bdf["len"].to_dict()

    # all backbone and edges
    df = pd.read_csv(args.edge_len_df, index_col=0)

    # find edges positions
    edges_pos = {}
    for e_id in df.columns.to_list():
        e = ut.Edge.from_str_id(e_id)
        # list of isolates that have this edge
        isolates = df[e_id].dropna().index.to_list()
        edges_pos[e_id] = find_edges_pos(pan, e, isolates)

    # write to file
    with open(args.positions, "w") as f:
        json.dump(edges_pos, f)
