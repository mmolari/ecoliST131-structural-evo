import numpy as np
import pandas as pd
import argparse
import pypangraph as ppg
from pypangraph.pangraph_projector import PanProjector
from itertools import combinations
from collections import defaultdict


# class pan_link:
#     """Class that represent an orientation-invariant link between two blocks,
#     that accounts for strandedness. If a and b are two links with the same topology:
#     - a == b
#     - hash(a) == hash(b)
#     - a in [b]
#     """

#     def __init__(self, bl_a, str_a, bl_b, str_b):
#         self.ida = bl_a
#         self.sa = str_a
#         self.idb = bl_b
#         self.sb = str_b

#     def __hash__(self) -> int:
#         hf = hash((self.ida, self.sa, self.idb, self.sb))
#         hr = hash((self.idb, not self.sb, self.ida, not self.sa))
#         return hf ^ hr

#     def __eq__(self, __o: object) -> bool:
#         res = (self.ida, self.sa, self.idb, self.sb) == (
#             __o.ida,
#             __o.sa,
#             __o.idb,
#             __o.sb,
#         )
#         res |= (self.idb, not self.sb, self.ida, not self.sa) == (
#             __o.ida,
#             __o.sa,
#             __o.idb,
#             __o.sb,
#         )
#         return res


def parse_args():
    parser = argparse.ArgumentParser(
        description="""extract pairwise distances from pangenome graph."""
    )
    parser.add_argument("--pangraph", type=str, help="input pangraph file")
    parser.add_argument("--csv", type=str, help="output distance dataframe")
    return parser.parse_args()


def pangraph_to_dist_df(pangraph_file):
    """Given a pangenome graph file, returns a dataframe of distances for all
    pairwise projections."""

    # load pangenome graph
    pan = ppg.Pangraph.load_json(pangraph_file)

    # set of strains
    strains = pan.strains()

    # create projector
    ppj = PanProjector(pan)

    df = []
    # cycle over all pairs of strains
    for acc_i, acc_j in combinations(strains, 2):

        # project over the pair
        pr = ppj.project(acc_i, acc_j)

        # n. blocks in projection
        n_blocks = len(set(pr.MPA.chunk_id).union(set(pr.MPB.chunk_id)))

        # block sequence length dictionary (tot alignment sequence in block)
        block_len = defaultdict(int)

        L_priv = 0
        L_shared = 0
        n_breakpoints = 0
        # for each of the two strains
        for M in [pr.MPA, pr.MPB]:

            s = M.comm  # mask for shared blocks
            n_breakpoints += np.sum(s != np.roll(s, 1))  # count number of breakpoints

            L = M.bl_Ls  # list of block lengths
            L_priv += np.sum(L[~s])
            L_shared += np.sum(L[s])

            cid = M.chunk_id
            for c, l in zip(cid, L):
                block_len[c] += l

        # evaluate partition entropy
        L_genomes = sum(block_len.values())
        entropy = (
            -sum(l * np.log(l / L_genomes) for l in block_len.values()) / L_genomes
        )

        # create dataframe entry
        for i, j in [[acc_i, acc_j], [acc_j, acc_i]]:
            res = {
                "si": i,
                "sj": j,
                "private seq. (bp)": L_priv,
                "shared seq. (bp)": L_shared,
                "n. breakpoints": n_breakpoints,
                "part. entropy": entropy,
                "n. blocks": n_blocks,
            }
            df.append(res)

    # add single strains
    for acc_i in strains:
        res = {
            "si": acc_i,
            "sj": acc_i,
            "private seq. (bp)": 0,
            "n. breakpoints": 0,
            "part. entropy": 0.0,
            "n. blocks": 1,
        }
        df.append(res)

    df = pd.DataFrame(df).sort_values(["si", "sj"])
    return df


if __name__ == "__main__":

    args = parse_args()

    # pangraph pairwise distances
    df = pangraph_to_dist_df(args.pangraph)

    # save dataframe, adding information header
    df_txt = df.to_csv(args.csv, index=False)
