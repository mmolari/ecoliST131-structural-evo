import argparse
import pypangraph as pp
import numpy as np
import pandas as pd

import utils as ut

from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Given a pangraph, generates a list of core edges and their frequencies.
        Only core-blocks longer than the specified threshold are considered.
        """
    )
    parser.add_argument("--pangraph", type=str, required=True)
    parser.add_argument("--csv", type=str, required=True)
    parser.add_argument("--len_thr", type=int, required=True)
    return parser.parse_args()


def edge_to_str(edge):
    "string representation for and edge, with alphabetical order"
    if edge.left.id < edge.right.id:
        return edge.to_str_id()
    else:
        return edge.invert().to_str_id()


if __name__ == "__main__":
    args = parse_args()

    # load pangraph
    pan = pp.Pangraph.load_json(args.pangraph)
    bdf = pan.to_blockstats_df()
    is_core = bdf["core"].to_dict()
    block_len = bdf["len"].to_dict()

    # turn to paths
    paths = ut.pangraph_to_path_dict(pan)

    # filter paths
    def core_filter(bid):
        if not is_core[bid]:
            return False
        if block_len[bid] < args.len_thr:
            return False
        return True

    core_paths = ut.filter_paths(paths, core_filter)

    # core-edges counter
    Edges = []
    for iso, path in core_paths.items():
        next_path = np.roll(path.nodes, -1)
        for l, r in zip(path.nodes, next_path):
            edge = ut.Edge(l, r)
            Edges.append(edge)
    Edges = Counter(Edges)

    # create and save dataframe
    df = pd.DataFrame(Edges.most_common(), columns=["edge", "count"])
    df["edge"] = df["edge"].apply(edge_to_str)
    df.to_csv(args.csv, index=False)
