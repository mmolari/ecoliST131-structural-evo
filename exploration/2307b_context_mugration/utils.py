import pypangraph as pp
import pathlib
from Bio import Phylo

LEN_THR = 500

root = pathlib.Path(__file__).parent.parent.parent.absolute()
pg_fld = root / "results" / "ST131" / "pangraph"
pan_file = pg_fld / "asm20-100-5-polished.json"
tree_file = pg_fld / "asm20-100-5-filtered-coretree.nwk"

expl_fld = root / "exploration" / "2307b_context_mugration" / "data"
named_nodes_tree_file = expl_fld / "named_tree.nwk"

fig_fld = root / "exploration" / "2307b_context_mugration" / "figs"


def load_pangraph():
    return pp.Pangraph.load_json(pan_file)


def load_tree():
    return Phylo.read(named_nodes_tree_file, "newick")


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

    def to_tuple(self):
        return (self.id, self.strand)

    @staticmethod
    def from_tuple(t) -> "Node":
        return Node(t[0], t[1])


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

    def to_list(self):
        return [n.to_tuple() for n in self.nodes]

    @staticmethod
    def from_list(path_list) -> "Path":
        return Path([Node.from_tuple(t) for t in path_list])


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


class Junction:
    """A junction is a combination of a node (or path) flanked by two nodes, with reverse-complement simmetry"""

    def __init__(self, left: Node, center, right: Node) -> None:
        self.left = left
        self.center = center
        self.right = right

    def invert(self) -> "Junction":
        return Junction(self.right.invert(), self.center.invert(), self.left.invert())

    def flanks_bid(self, bid) -> bool:
        return (self.left.id == bid) or (self.right.id == bid)

    def __side_eq__(self, o: object) -> bool:
        return self.left == o.left and self.center == o.center and self.right == o.right

    def __eq__(self, o: object) -> bool:
        return self.__side_eq__(o) or self.__side_eq__(o.invert())

    def __side_hash__(self) -> int:
        return hash((self.left, self.center, self.right))

    def __hash__(self) -> int:
        return self.__side_hash__() ^ self.invert().__side_hash__()

    def __repr__(self) -> str:
        return f"{self.left} <-- {self.center} --> {self.right}"


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


def to_core_adjacencies(nodes, is_core):
    """Given a list of nodes, links each accessory block to the flanking core blocks."""
    adj = []
    N = len(nodes)

    for i, n in enumerate(nodes):
        if is_core[n.id]:
            continue
        bef, aft = None, None
        bi, ai = i, i
        while bef is None:
            bi = (bi - 1) % N
            if is_core[nodes[bi].id]:
                bef = nodes[bi]
        while aft is None:
            ai = (ai + 1) % N
            if is_core[nodes[ai].id]:
                aft = nodes[ai]
        j = Junction(bef, n, aft)
        adj.append(j)

    return adj
