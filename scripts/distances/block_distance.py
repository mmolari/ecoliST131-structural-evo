import argparse

import numpy as np
import pandas as pd
import pypangraph as pp


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
    parser.add_argument("--len_thr", help="block length threshold", type=int)
    return parser.parse_args()


if __name__ == "__main__":
    # load pangraph
    args = parse_args()
    pan = pp.Pangraph.load_json(args.pan)

    strains = np.sort(pan.strains())
    N = len(strains)

    # block count matrix to PA
    df = (pan.to_blockcount_df() > 0).astype(int)

    propr = pan.to_blockstats_df()
    mask = propr["len"] > args.len_thr
    df = df.loc[:, mask]

    mask &= ~propr["core"] & ~propr["duplicated"]
    acc_blocks = propr[mask].index
    acc_df = df.loc[:, acc_blocks]

    # evaluate distance and sharing matrices
    M_pa, M_s, M_pa_acc = [np.zeros((N, N), int) for _ in range(3)]
    for i in range(N):
        si = strains[i]
        pi = df.loc[si]
        M_s[i, i] = pi.sum()
        for j in range(i + 1, N):
            sj = strains[j]
            pj = df.loc[sj]
            n_shared = (pi & pj).sum()
            n_diff = (pi ^ pj).sum()
            n_diff_acc = (acc_df.loc[si] ^ acc_df.loc[sj]).sum()

            M_s[i, j] = n_shared
            M_s[j, i] = n_shared
            M_pa[i, j] = n_diff
            M_pa[j, i] = n_diff
            M_pa_acc[i, j] = n_diff_acc
            M_pa_acc[j, i] = n_diff_acc

    # matrix to dataframe
    df_pa = matrix_to_df(M_pa, S=strains, sname="block_PA")
    df_s = matrix_to_df(M_s, S=strains, sname="block_sharing")
    df_pa_acc = matrix_to_df(M_pa_acc, S=strains, sname="acc_block_PA")

    # concatenate and save
    df = pd.concat([df_pa, df_s, df_pa_acc], axis=1, verify_integrity=True)
    df.to_csv(args.csv)
