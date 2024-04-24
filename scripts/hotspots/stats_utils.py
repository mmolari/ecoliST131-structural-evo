from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from copy import deepcopy
import itertools


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
        for i in range(L - 1):
            e = Edge(p.nodes[i], p.nodes[i + 1])
            ct.update([e])
    return dict(ct)


def path_block_count(paths):
    """Count internal blocks of paths"""
    ct = Counter()
    for iso, p in paths.items():
        for node in p.nodes:
            ct.update([node.id])
    return dict(ct)


def core_alignment(pan):
    """Returns a matrix of gapless core SNPs and total alignment length."""
    strains = pan.strains()
    bdf = pan.to_blockstats_df()
    core_blocks = bdf.query("core == True").index
    As = defaultdict(str)
    for bid in core_blocks:
        b = pan.blocks[bid]
        A, O = b.alignment.generate_alignments()
        for a, o in zip(A, O):
            As[o[0]] += a
    # turn into alignment matrix:
    As = np.array([np.array(list(As[s])) for s in strains])
    # remove gap columns
    mask = np.all(As != "-", axis=0)
    As = As[:, mask]
    return As


def print_progress_bar(k, N):
    end = "\n" if k == N - 1 else "\r"
    fr = (k + 1) / N
    print(f"[{'#' * int(50 * fr):50s}] {100 * fr:.2f}%", end=end)


def core_divergence(A, si, sj, strains):
    S = list(strains)
    i = S.index(si)
    j = S.index(sj)
    Ai = A[i, :]
    Aj = A[j, :]
    div = np.mean(Ai != Aj)
    return {"local_core_div": div, "alignment_len": len(Ai)}


def block_distance_stats(si, sj, bpa, bdf):
    bi = bpa.loc[si]
    bj = bpa.loc[sj]
    private_blocks = np.abs(bi - bj)
    private_blocks_len = (bdf.loc[private_blocks.index]["len"] * private_blocks).sum()
    shared_blocks = np.minimum(bi, bj)
    shared_blocks_len = (bdf.loc[shared_blocks.index]["len"] * shared_blocks).sum()
    return {
        "n. private blocks": private_blocks.sum(),
        "len. private blocks": private_blocks_len,
        "n. shared blocks": shared_blocks.sum(),
        "len. shared blocks": shared_blocks_len,
    }


def breakpoint_stats(si, sj, paths):
    subpaths = {iso: paths[iso] for iso in [si, sj]}

    edge_ct = path_edge_count(subpaths)
    block_ct = path_block_count(subpaths)

    n_breakpoints = 0
    for e, ec in edge_ct.items():
        bl, br = e.left.id, e.right.id
        if (ec == block_ct[bl]) and (ec == block_ct[br]):
            continue
        n_breakpoints += 1

    return {
        "n. breakpoints": n_breakpoints,
    }


class MergePath:
    def __init__(self, subpaths, lengths, is_core) -> None:
        self.subpaths = subpaths
        self.lengths = lengths
        self.is_core = is_core

    def __remove_from_path(self, edge, acc) -> None:
        for k, sp in self.subpaths.items():
            L = len(sp.nodes)
            to_remove = None
            for l in range(L - 1):
                l_edge = Edge(sp.nodes[l], sp.nodes[l + 1])
                if l_edge == edge:
                    if sp.nodes[l].id == acc:
                        to_remove = l
                    elif sp.nodes[l + 1].id == acc:
                        to_remove = l + 1
                    break
            assert to_remove is not None
            sp.nodes.pop(to_remove)
            assert len(sp.nodes) == L - 1

    def merge(self, edge, core, acc) -> None:
        self.__remove_from_path(edge, acc)
        self.lengths[core] += self.lengths[acc]

    def recursive_extend(self) -> None:
        changed = False
        edge_ct = path_edge_count(self.subpaths)
        for e, ec in edge_ct.items():
            bl, br = e.left.id, e.right.id
            cl, cr = self.is_core[bl], self.is_core[br]
            if cl and (ec == len(self.subpaths)):
                # print(f"merging {bl} <- {br}")
                self.merge(e, bl, br)
                changed = True
                break
            elif cr and (ec == len(self.subpaths)):
                # print(f"merging {br} <- {bl}")
                self.merge(e, br, bl)
                changed = True
                break
        if changed:
            self.recursive_extend()

    def stats(self) -> dict:
        n_edges = len(path_edge_count(self.subpaths))
        block_ct = path_block_count(self.subpaths)
        core_len, acc_len = 0, 0
        n_shared, n_private = 0, 0
        for bid, n in block_ct.items():
            if self.is_core[bid]:
                core_len += self.lengths[bid]
                n_shared += 1
            else:
                acc_len += self.lengths[bid]
                n_private += 1
        return {
            "merge_n_edges": n_edges,
            "merge_core_len": core_len,
            "merge_acc_len": acc_len,
            "merge_n_shared": n_shared,
            "merge_n_private": n_private,
        }


def merge_paths_stats(si, sj, paths, bdf):
    bl_lens = bdf["len"].to_dict()
    bl_core = bdf["core"].to_dict()
    subpaths = {iso: deepcopy(paths[iso]) for iso in [si, sj]}
    mp = MergePath(subpaths, bl_lens, bl_core)
    mp.recursive_extend()
    return mp.stats()


def extract_hotspot_stats(pan):

    df = []

    # calculate utilities
    bdf = pan.to_blockstats_df()
    bpa = pan.to_blockcount_df()
    strains = pan.strains()
    path_dict = pangraph_to_path_dict(pan)
    A = core_alignment(pan)

    # iterate over all strain pairs
    N = len(strains)
    for k, (si, sj) in enumerate(itertools.combinations(strains, 2)):

        if k % 100 == 0:
            print_progress_bar(k, N * (N - 1) // 2)

        # ensure si > sj
        if si < sj:
            si, sj = sj, si
        res = {"si": si, "sj": sj}

        res |= block_distance_stats(si, sj, bpa, bdf)
        res |= breakpoint_stats(si, sj, path_dict)
        res |= core_divergence(A, si, sj, strains)
        res |= merge_paths_stats(si, sj, path_dict, bdf)

        df.append(res)
    df = pd.DataFrame(df)
    df.set_index(["si", "sj"], inplace=True)
    return df
