# %%
from collections import defaultdict, Counter


class Node:
    """Combination of block id and strandedness"""

    def __init__(self, bid: str, strand: bool) -> None:
        self.id = bid
        self.strand = strand

    def invert(self) -> "Node":
        return Node(self.id, not self.strand)

    def __eq__(self, other: object) -> bool:
        return self.id == other.id and self.strand == other.strand

    def __hash__(self) -> int:
        return hash((self.id, self.strand))

    def __repr__(self) -> str:
        s = "+" if self.strand else "-"
        return f"[{self.id}|{s}]"

    def to_str_id(self):
        s = "f" if self.strand else "r"
        return f"{self.id}_{s}"

    @staticmethod
    def from_str_id(t) -> "Node":
        bid = t.split("_")[0]
        strand = True if t.split("_")[1] == "f" else False
        return Node(bid, strand)


class Path:
    """A path is a list of nodes"""

    def __init__(self, nodes=[]) -> None:
        self.nodes = nodes

    def add_left(self, node: Node) -> None:
        self.nodes.insert(0, node)

    def add_right(self, node: Node) -> None:
        self.nodes.append(node)

    def invert(self) -> "Path":
        return Path([n.invert() for n in self.nodes[::-1]])

    def __eq__(self, o: object) -> bool:
        return self.nodes == o.nodes

    def __hash__(self) -> int:
        return hash(tuple(self.nodes))

    def __repr__(self) -> str:
        return "_".join([str(n) for n in self.nodes])

    def __len__(self) -> int:
        return len(self.nodes)

    def to_list(self):
        return [n.to_str_id() for n in self.nodes]

    @staticmethod
    def from_list(path_list) -> "Path":
        return Path([Node.from_str_id(nid) for nid in path_list])


class Edge:
    """Oriented link between two nodes/paths"""

    def __init__(self, left, right) -> None:
        self.left = left
        self.right = right

    def invert(self) -> "Edge":
        return Edge(self.right.invert(), self.left.invert())

    def __side_eq__(self, o: object) -> bool:
        return self.left == o.left and self.right == o.right

    def __eq__(self, o: object) -> bool:
        return self.__side_eq__(o) or self.__side_eq__(o.invert())

    def __side_hash__(self) -> int:
        return hash((self.left, self.right))

    def __hash__(self) -> int:
        return self.__side_hash__() ^ self.invert().__side_hash__()

    def __repr__(self) -> str:
        return f"{self.left} <--> {self.right}"

    def __to_str_id(self) -> str:
        return "__".join([self.left.to_str_id(), self.right.to_str_id()])

    def to_str_id(self) -> str:
        A = self.__to_str_id()
        B = self.invert().__to_str_id()
        return A if A < B else B

    @staticmethod
    def from_str_id(t) -> "Edge":
        left, right = t.split("__")
        return Edge(Node.from_str_id(left), Node.from_str_id(right))


def pangraph_to_path_dict(pan):
    """Creates a dictionary isolate -> path objects"""
    res = {}
    for path in pan.paths:
        name = path.name
        B = path.block_ids
        S = path.block_strands
        nodes = [Node(b, s) for b, s in zip(B, S)]
        res[name] = Path(nodes)
    return res


def filter_paths(paths, keep_f):
    """Given a filter function, removes nodes that fail the condition from
    the path dictionaries."""
    res = {}
    for iso, path in paths.items():
        filt_path = Path([node for node in path.nodes if keep_f(node.id)])
        res[iso] = filt_path
    return res


def path_categories(paths):
    """Returns a list of touples, one per non-empty path, with the following info:
    (count, path, [list of isolates])"""
    iso_list = defaultdict(list)
    n_paths = defaultdict(int)
    nodes = {}
    for iso, path in paths.items():
        if len(path.nodes) > 0:
            n_paths[path] += 1
            iso_list[path].append(iso)
            nodes[path] = path.nodes

    # sort by count
    path_cat = [(count, nodes[path], iso_list[path]) for path, count in n_paths.items()]
    path_cat.sort(key=lambda x: x[0], reverse=True)
    return path_cat


def path_edge_count(paths):
    """Count internal edges of paths"""
    ct = Counter()
    for iso, p in paths.items():
        L = len(p.nodes)
        for i in range(L):
            e = Edge(p.nodes[i], p.nodes[(i + 1) % L])
            ct.update([e])
    return dict(ct)


# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import pathlib
import pypangraph as pp
from Bio import SeqIO

dest = pathlib.Path("data/single_block")
dest.mkdir(exist_ok=True, parents=True)

# %%
fname = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(fname)
bdf = pan.to_blockstats_df()
# %%
paths = pangraph_to_path_dict(pan)
paths = filter_paths(paths, lambda bid: bdf.loc[bid, "core"])
# %%
ref = "NZ_CP096110.1"
qry = "NZ_CP059281.1"
qry = "NZ_CP059279.1"


pr = paths[ref]
pq = paths[qry]

ec = path_edge_count({"r": pr, "q": pq})
ec = {k: v for k, v in ec.items() if v == 1}
ns = [e.left.id for e in ec.keys()] + [e.right.id for e in ec.keys()]
ns = Counter(ns)


# %%
probl_block = "HEVMPTBCLE"
seq = pan.blocks[probl_block].sequence
seq = SeqIO.SeqRecord(seq, id=probl_block, description="")
bfile = dest / f"{probl_block}.fa"
SeqIO.write([seq], bfile, "fasta")


# %%
ref_fname = f"../../data/fa/{ref}.fa"
qry_fname = f"../../data/fa/{qry}.fa"

# map with minimap2
out_fname = dest / f"{probl_block}_{qry}.paf"
cmd = f"minimap2 -x asm5 {bfile} {qry_fname} > {out_fname}"
subprocess.run(cmd, shell=True)

out_fname = dest / f"{probl_block}_{ref}.paf"
cmd = f"minimap2 -x asm5 {bfile} {ref_fname} > {out_fname}"
subprocess.run(cmd, shell=True)

# %%
