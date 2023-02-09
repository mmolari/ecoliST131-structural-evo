import treetime
import argparse
import json
from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser(
        description="refine core genome tree using treetime"
    )
    parser.add_argument("--tree_in", help="input nwk file", type=str)
    parser.add_argument("--json", help="reduced alignment information", type=str)
    parser.add_argument("--aln", help="alignment fasta file", type=str)
    parser.add_argument("--tree_out", help="output nwk file", type=str)
    args = parser.parse_args()
    return args


def rescale_branch_length(node, factor):
    if node.branch_length is not None:
        node.branch_length *= factor
    for child in node:
        rescale_branch_length(child, factor)


if __name__ == "__main__":

    args = parse_args()

    # load alignment information and evaluate rescaling factor
    with open(args.json, "r") as f:
        info = json.load(f)
    aln_L = info["n. consensus"] + info["n. snps"]
    factor = info["n. snps"] / aln_L

    # load input tree, midpoint root and rescale
    tree = Phylo.read(args.tree_in, format="newick")
    rescale_branch_length(tree.root, factor)
    tree.root_at_midpoint()
    tree.ladderize()

    # instantiate treetime
    myTree = treetime.TreeAnc(
        gtr="Jukes-Cantor",
        tree=tree,
        aln=args.aln,
        verbose=0,
        seq_len=aln_L,
    )

    myTree.tree.root.branch_length = 0.0
    # optimize branch length
    print("--------- Treetime ---------")
    print("    < before optimizing >")
    print("tot. branch length:", myTree.tree.total_branch_length())
    print("n. nonterminals:", len(myTree.tree.get_nonterminals()))
    myTree.optimize_tree(prune_short=True)
    print("    < after optimizing >")
    print("tot. branch length:", myTree.tree.total_branch_length())
    print("n. nonterminals:", len(myTree.tree.get_nonterminals()))

    # save tree
    Phylo.write(
        myTree.tree,
        args.tree_out,
        format="newick",
        # format_branch_length="%1.10f",
        format_branch_length="%.5e",
    )
