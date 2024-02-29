import pandas as pd
import pypangraph as pp
import utils as ut
from collections import defaultdict
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="""Giv.""")
    parser.add_argument("--pangraph", type=str, required=True)
    parser.add_argument("--edge_df", type=str, required=True)
    parser.add_argument("--out_edge_len_df", type=str, required=True)
    return parser.parse_args()


def load_pangraph_data(pangraph_fname):
    pan = pp.Pangraph.load_json(pangraph_fname)
    path_dict = ut.pangraph_to_path_dict(pan)
    block_L = pan.to_blockstats_df()["len"].to_dict()
    return path_dict, block_L


def load_edge_data(edge_df_fname):
    edf = pd.read_csv(edge_df_fname, index_col=0)
    edf = edf > 0
    edges = [ut.Edge.from_str_id(e) for e in edf.columns]
    core_blocks = [e.left.id for e in edges] + [e.right.id for e in edges]
    core_blocks = list(set(core_blocks))
    return edf, edges, core_blocks


def extract_edge_acc_blocks(path_dict, edges, edf):
    def is_coreblock(node_id):
        return node_id in core_blocks

    edge_acc_blocks = defaultdict(set)
    for iso, path in path_dict.items():
        junctions = ut.path_junction_split(path, is_coreblock)
        for J in junctions:
            edge = J.flanking_edge()
            assert edge in edges, f"{edge} not in {edges}"
            if len(J.center.nodes) == 0:
                continue
            assert edf.loc[iso][edge.to_str_id()], f"{edge} not in {iso}"
            edge_blocks = [n.id for n in J.center.nodes]
            edge_acc_blocks[edge].update(edge_blocks)
    return edge_acc_blocks


if __name__ == "__main__":

    args = parse_args()

    path_dict, block_L = load_pangraph_data(args.pangraph)
    edf, edges, core_blocks = load_edge_data(args.edge_df)

    edge_acc_blocks = extract_edge_acc_blocks(path_dict, edges, edf)

    edge_acc_pangenome = {
        e.to_str_id(): (sum(block_L[b] for b in blocks), len(blocks))
        for e, blocks in edge_acc_blocks.items()
    }

    edf_len = pd.DataFrame(edge_acc_pangenome).T
    edf_len.columns = ["pangnome_len", "pangenome_n_blocks"]

    edf_len.to_csv(args.out_edge_len_df)
