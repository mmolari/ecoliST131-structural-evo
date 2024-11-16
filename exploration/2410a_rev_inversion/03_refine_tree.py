import treetime
import argparse
import pandas as pd
import numpy as np
from Bio import Phylo, AlignIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="refine core genome tree using treetime"
    )
    parser.add_argument("--tree_in", help="input nwk file", type=str)
    parser.add_argument("--aln", help="reduced alignment fasta file", type=str)
    parser.add_argument("--lengths", help="csv file with alignment lengths", type=str)
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

    # load input tree, midpoint root and rescale
    tree = Phylo.read(args.tree_in, format="newick")
    tree.root_at_midpoint()
    tree.ladderize()

    # load lengths
    ldf = pd.read_csv(args.lengths)
    full_len = ldf["L_ungapped"][0]
    restr_len = ldf["L_restr"][0]
    # L_factor = restr_len / full_len

    # instantiate treetime
    myTree = treetime.TreeAnc(
        gtr="Jukes-Cantor",
        tree=tree,
        aln=args.aln,
        verbose=0,
        seq_len=full_len,
    )

    myTree.tree.root.branch_length = 0.0
    # optimize branch length
    print("--------- Treetime ---------")
    print("    < before optimizing >")
    print("tot. branch length:", myTree.tree.total_branch_length())
    print("n. nonterminals:", len(myTree.tree.get_nonterminals()))
    myTree.optimize_tree()
    # rescale_branch_length(myTree.tree.root, L_factor)
    print("    < after optimizing >")
    print("tot. branch length:", myTree.tree.total_branch_length())
    print("n. nonterminals:", len(myTree.tree.get_nonterminals()))

    # save tree
    myTree.tree.root_at_midpoint()
    myTree.tree.ladderize()
    Phylo.write(
        myTree.tree,
        args.tree_out,
        format="newick",
        # format_branch_length="%1.10f",
        format_branch_length="%.5e",
    )
