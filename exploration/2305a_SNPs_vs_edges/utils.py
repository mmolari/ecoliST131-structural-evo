import numpy as np
from Bio import AlignIO, Phylo
import pypangraph as pp
import pandas as pd
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


def path_to_blocklist(pan, k, mask_dict):
    p = pan.paths[k]
    return [
        Block(bid, s) for bid, s in zip(p.block_ids, p.block_strands) if mask_dict[bid]
    ]


def blocklist_to_edgelist(bl):
    return [Edge(b1, b2) for b1, b2 in zip(bl, np.roll(bl, -1))]


def extract_edge_dictionary(pangraph_file, len_thr=0, remove_dupl=False):
    pan = pp.Pangraph.load_json(pangraph_file)

    propr = pan.to_blockstats_df()
    block_mask = propr["len"] > len_thr
    if remove_dupl:
        block_mask &= ~propr["duplicated"]
    # block_mask &= propr["core"]
    block_mask_dict = block_mask.to_dict()

    # extract path and edge representation
    strains = set(pan.strains())
    paths = {k: path_to_blocklist(pan, k, block_mask_dict) for k in strains}
    edges = {k: blocklist_to_edgelist(bl) for k, bl in paths.items()}

    #  dictionary {edge -> strain set}
    edge_strains = defaultdict(set)
    for k, el in edges.items():
        for e in el:
            edge_strains[e].add(k)

    return edge_strains


def extract_block_PA_df(pangraph_file, len_thr=0, remove_dupl=False):
    pan = pp.Pangraph.load_json(pangraph_file)

    # block count matrix to PA
    df = (pan.to_blockcount_df() > 0).astype(int)
    propr = pan.to_blockstats_df()
    mask = propr["len"] > len_thr
    if remove_dupl:
        mask = mask & (~propr["duplicated"])
    df = df.loc[:, mask]

    return df


def aln_to_branchings(fasta_file):
    aln = AlignIO.read(fasta_file, "fasta")
    A = np.array(aln)
    strains = [s.name for s in aln]

    N, L = A.shape

    branchings = []
    for l in range(L):
        row = A[:, l]

        # check biallelic
        alleles, cts = np.unique(row, return_counts=True)

        if len(alleles) > 2:
            print(f"Site {l} has more than 2 alleles")
            continue

        if len(alleles) < 2:
            print(f"Site {l} has less than 2 alleles")
            continue

        # find branching
        allele = alleles[np.argmin(cts)]
        idxs = np.argwhere(row == allele).flatten()

        # save pair
        branchings.append([strains[i] for i in idxs])

    return branchings
