import json
import argparse

import pypangraph as pp
import pandas as pd

import utils as ut


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Given a pangraph, generates a list of core edges and their frequencies.
        Only core-blocks longer than the specified threshold are considered.
        """
    )
    parser.add_argument("--pangraph", type=str, required=True)
    parser.add_argument("--core_edges_csv", type=str, required=True)
    parser.add_argument("--paths_json", type=str, required=True)
    return parser.parse_args()


def to_core_junctions(path, backbone_ids):
    """Given a path, returns a dictionary of backbone edges -> junctions"""
    left_core = None
    right_core = None
    first_core = None
    CJs = {}
    center = []

    # set index to first backbone node
    i = 0
    while True:
        node = path.nodes[i]
        i += 1
        is_backbone = node.id in backbone_ids
        if is_backbone:
            first_core = node
            left_core = node
            break

    L = len(path.nodes)
    while True:
        node = path.nodes[i]
        is_backbone = node.id in backbone_ids
        if is_backbone:
            right_core = node
            J = ut.Junction(left_core, ut.Path(center), right_core)
            E = ut.Edge(left_core, right_core)
            CJs[E] = J
            center = []
            if right_core == first_core:
                break
            else:
                left_core = right_core
        else:
            center.append(node)

        i = (i + 1) % L

    return CJs


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
    df = pd.read_csv(args.core_edges_csv)
    df = df[df["count"] == N_iso]
    backbone_edges = [ut.Edge.from_str_id(e) for e in df["edge"]]
    backbone_blocks = set(
        [e.left.id for e in backbone_edges] + [e.right.id for e in backbone_edges]
    )

    # extract paths
    paths = ut.pangraph_to_path_dict(pan)

    # extract backbone junctions
    BBj = {}
    for iso, path in paths.items():
        BBj[iso] = to_core_junctions(path, backbone_blocks)

    # swap keys of nested dictionary: first edges and then isolates
    backbone_junctions = {}
    for e in backbone_edges:
        key = e.to_str_id()
        backbone_junctions[key] = {}
        for iso in strains:
            backbone_junctions[key][iso] = BBj[iso][e].to_list()

    # save in json file
    with open(args.paths_json, "w") as f:
        json.dump(backbone_junctions, f)
