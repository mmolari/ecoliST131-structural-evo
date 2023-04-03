import copy
import argparse

import numpy as np
import pandas as pd
import pypangraph as pp

from collections import defaultdict


def matrix_to_df(M, S, sname):
    si = pd.Series(S, name="si")
    sj = pd.Series(S, name="sj")
    df = pd.DataFrame(M, index=si, columns=sj)
    series = df.unstack()
    series.name = sname
    return pd.DataFrame(series)


def parse_args():
    parser = argparse.ArgumentParser(description="Evalaute edge-related distances")
    parser.add_argument("--pan", help="pangenome graph", type=str)
    parser.add_argument("--csv", help="output dataframe", type=str)
    return parser.parse_args()


if __name__ == "__main__":

    # load pangraph
    args = parse_args()
    pan = pp.Pangraph.load_json(args.pan)

    strains = np.sort(pan.strains())
    N = len(strains)

    # block count matrix to PA
    df = (pan.to_blockcount_df() > 0).astype(int)

    # evaluate distance and sharing matrices
    M_pa, M_s = [np.zeros((N, N), int) for _ in range(2)]
    for i in range(N):
        si = strains[i]
        pi = df.loc[si]
        M_s[i, i] = pi.sum()
        for j in range(i + 1, N):
            sj = strains[j]
            pj = df.loc[sj]
            n_shared = (pi & pj).sum()
            n_diff = (pi ^ pj).sum()

            M_s[i, j] = n_shared
            M_s[j, i] = n_shared
            M_pa[i, j] = n_diff
            M_pa[j, i] = n_diff

    # matrix to dataframe
    df_pa = matrix_to_df(M_pa, S=strains, sname="block_PA")
    df_s = matrix_to_df(M_s, S=strains, sname="block_sharing")

    # concatenate and save
    df = pd.concat([df_pa, df_s], axis=1, verify_integrity=True)
    df.to_csv(args.csv)
