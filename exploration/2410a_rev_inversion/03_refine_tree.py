import treetime
import argparse
import numpy as np
from Bio import Phylo, AlignIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="refine core genome tree using treetime"
    )
    parser.add_argument("--tree_in", help="input nwk file", type=str)
    parser.add_argument("--aln", help="reduced alignment fasta file", type=str)
    parser.add_argument("--tree_out", help="output nwk file", type=str)
    args = parser.parse_args()
    return args


def aln_len(aln_file):
    aln = AlignIO.read(aln_file, format="fasta")
    A = np.array(aln)
    assert np.all(A != "-")
    return A.shape[1]


if __name__ == "__main__":
    args = parse_args()

    tot_len = aln_len(args.aln)

    # load input tree, midpoint root and rescale
    tree = Phylo.read(args.tree_in, format="newick")
    tree.root_at_midpoint()
    tree.ladderize()

    # instantiate treetime
    myTree = treetime.TreeAnc(
        gtr="Jukes-Cantor",
        tree=tree,
        aln=args.aln,
        verbose=0,
        seq_len=tot_len,
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
    myTree.tree.root_at_midpoint()
    myTree.tree.ladderize()
    Phylo.write(
        myTree.tree,
        args.tree_out,
        format="newick",
        # format_branch_length="%1.10f",
        format_branch_length="%.5e",
    )
