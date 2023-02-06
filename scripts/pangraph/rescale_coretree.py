import json
import argparse
from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser(description="rescale core genome tree")
    parser.add_argument("--tree_in", help="input nwk file", type=str)
    parser.add_argument("--tree_out", help="output nwk file", type=str)
    parser.add_argument("--json", help="reduced alignment information", type=str)
    parser.add_argument(
        "--confidence",
        help="collapses branches with below-threshold confidence",
        type=float,
    )
    args = parser.parse_args()
    return args


def rescale_branch_length(node, factor):
    if node.branch_length is not None:
        node.branch_length *= factor
    for child in node:
        rescale_branch_length(child, factor)


if __name__ == "__main__":

    args = parse_args()

    # evaluate rescale factor
    with open(args.json, "r") as f:
        info = json.load(f)
    factor = info["n. snps"] / info["n. consensus"]

    # load tree
    tree = Phylo.read(args.tree_in, format="newick")
    # collapse below-confidence nodes
    tree.collapse_all(
        lambda c: c.confidence is not None and c.confidence < args.confidence
    )

    # root at midpoint and ladderize
    tree.root_at_midpoint()
    tree.ladderize()

    # rescale
    rescale_branch_length(tree.root, factor)

    # write
    Phylo.write(tree, args.tree_out, format="newick")
