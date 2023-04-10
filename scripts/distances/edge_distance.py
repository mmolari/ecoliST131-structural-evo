import copy
import argparse

import numpy as np
import pandas as pd
import pypangraph as pp

from collections import defaultdict


class Block:
    def __init__(self, bid: str, s: bool):
        self.id = bid
        self.s = s

    def __repr__(self):
        sign = "+" if self.s else "-"
        return f"[{self.id}|{sign}]"

    def __eq__(self, b):
        return (self.id == b.id) and (self.s == b.s)

    def __hash__(self):
        return hash((self.id, self.s))

    def invert(self):
        return Block(self.id, not self.s)


class Edge:
    def __init__(self, b1: Block, b2: Block):
        self.b = (b1, b2)

    def __eq__(self, e):
        if self.b == e.b:
            return True
        elif self.invert().b == e.b:
            return True
        else:
            return False

    def __repr__(self):
        return f"({self.b[0]} - {self.b[1]})"

    def __hash__(self):
        e = self.invert()
        return hash((self.b[0], self.b[1])) ^ hash((e.b[0], e.b[1]))

    def invert(self):
        return Edge(self.b[1].invert(), self.b[0].invert())

    def has_block_from(self, bl_ids):
        return (self.b[0].id in bl_ids) or (self.b[1].id in bl_ids)


def path_to_blocklist(pan, k):
    p = pan.paths[k]
    return [Block(bid, s) for bid, s in zip(p.block_ids, p.block_strands)]


def blocklist_to_edgelist(bl):
    return [Edge(b1, b2) for b1, b2 in zip(bl, np.roll(bl, -1))]


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

    args = parse_args()

    pan = pp.Pangraph.load_json(args.pan)

    # extract path and edge representation
    strains = set(pan.strains())
    N = len(strains)
    Bids = pan.to_paths_dict()
    paths = {k: path_to_blocklist(pan, k) for k in strains}
    edges = {k: blocklist_to_edgelist(bl) for k, bl in paths.items()}

    #  dictionary {edge -> strain set}
    edge_strains = defaultdict(set)
    for k, el in edges.items():
        for e in el:
            edge_strains[e].add(k)

    N = len(strains)
    strain_order = sorted(list(strains))
    strain_pos = {s: n for n, s in enumerate(strain_order)}

    # edge P/A distance matrix and edge sharing matrix
    M_pa, M_s, M_pa_r = [np.zeros((N, N), int) for _ in range(3)]
    for e, s in edge_strains.items():
        for s1 in s:
            i1 = strain_pos[s1]
            for s2 in s:
                i2 = strain_pos[s2]
                M_s[i1, i2] += 1
            for s2 in strains - s:
                i2 = strain_pos[s2]
                M_pa[i1, i2] += 1
                M_pa[i2, i1] += 1
                if e.has_block_from(Bids[s2]):
                    M_pa_r[i1, i2] += 1
                    M_pa_r[i2, i1] += 1
    # matrix to dataframe
    df_pa = matrix_to_df(M_pa, S=strain_order, sname="edge_PA")
    df_pa_r = matrix_to_df(M_pa_r, S=strain_order, sname="edge_PA_reduced")
    df_s = matrix_to_df(M_s, S=strain_order, sname="edge_sharing")

    # concatenate and save
    df = pd.concat([df_pa, df_pa_r, df_s], axis=1, verify_integrity=True)
    df.to_csv(args.csv)
