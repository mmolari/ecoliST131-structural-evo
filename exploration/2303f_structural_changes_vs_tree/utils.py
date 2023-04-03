import numpy as np


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


def path_to_blocklist(pan, k):
    p = pan.paths[k]
    return [Block(bid, s) for bid, s in zip(p.block_ids, p.block_strands)]


def blocklist_to_edgelist(bl):
    return [Edge(b1, b2) for b1, b2 in zip(bl, np.roll(bl, -1))]


def dupl_blocks(pan):
    df = pan.to_blockstats_df()
    return df[df["duplicated"]].index.to_list()


def pangraph_to_block_splits(pan):
    df = pan.to_blockcount_df() > 0
    splits = {}
    for bl, pa in df.T.iterrows():
        splits[bl] = set(pa[pa].index.to_list())
    return splits
