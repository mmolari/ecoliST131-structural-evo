import argparse
import json
import numpy as np
import pandas as pd
from Bio import AlignIO

from itertools import combinations


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate pairwise distances from core alignment."
    )
    parser.add_argument("--aln", help="restricted fasta core alignment", type=str)
    parser.add_argument("--n_cons", help="n. consensus sites", type=int)
    parser.add_argument("--csv", help="output csv distance file", type=str)
    parser.add_argument("--label", help="distance prefix", type=str)
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # load alignment as matrix
    aln = AlignIO.read(args.aln, format="fasta")
    M = np.array(aln)

    # alignment length
    L = aln.get_alignment_length()

    #  evaluate divergence factor
    n_snps = L
    factor = n_snps / (n_snps + args.n_cons)

    # label
    lab = args.label

    # get strain names
    strains = [row.name for row in aln]

    # create distance dataframe
    df = []
    for si, sj in combinations(strains, 2):
        ni = strains.index(si)
        nj = strains.index(sj)
        d = np.sum(M[ni] != M[nj]) / L
        d *= factor
        df.append({"si": si, "sj": sj, lab: d})
        df.append({"si": sj, "sj": si, lab: d})

    for s in strains:
        df.append({"si": s, "sj": s, lab: 0.0})

    df = pd.DataFrame(df).sort_values(["si", "sj"])

    # save as csv
    df.to_csv(args.csv, index=False)

    # convert to distance matrix
    # D_mat = df.pivot(index="si", columns="sj", values="core_div").fillna(0)
