import pypangraph as pp
from Bio import Phylo
import pathlib
import unittest
import json
import itertools as itt


root = pathlib.Path(__file__).parent.parent.parent.absolute()
pg_fld = root / "results" / "ST131" / "pangraph"
expl_fld = root / "exploration" / "2307a_context_mugration" / "data"
pan_file = pg_fld / "asm20-100-5-polished.json"
tree_file = pg_fld / "asm20-100-5-filtered-coretree.nwk"
named_nodes_tree_file = expl_fld / "named_tree.nwk"


def load_pangraph():
    return pp.Pangraph.load_json(pan_file)


def load_tree():
    return Phylo.read(tree_file, "newick")


def load_nodenamed_tree():
    return Phylo.read(named_nodes_tree_file, "newick")


def to_core_adjacencies(node_path, is_core):
    adj = {}
    N = len(node_path)

    for i, n in enumerate(node_path):
        if is_core[n.id]:
            continue
        bef, aft = None, None
        bi, ai = i, i
        while bef is None:
            bi = (bi - 1) % N
            if is_core[node_path[bi].id]:
                bef = node_path[bi]
        while aft is None:
            ai = (ai + 1) % N
            if is_core[node_path[ai].id]:
                aft = node_path[ai]
        j = Junction(bef, n, aft)
        if not n.strand:
            j = j.invert()
        adj[n.id] = j

    return adj


def parse_gtr(fname):
    with open(fname, "r") as f:
        lines = f.readlines()

    info = {}

    def find_subsequent(beg):
        res = []
        keep = False
        for l in lines:
            if l.startswith(beg):
                keep = True
            if len(l.strip()) == 0:
                keep = False
            if keep:
                res.append(l.strip())
        return res

    def parse_matrix(mp, mat):
        c1 = mat[0].split("\t")
        c1 = [mp[c] for c in c1]
        c2 = [mp[m.split("\t")[0]] for m in mat[1:]]
        r1 = [float(m) for m in mat[1].split("\t")[1:]]
        r2 = [float(m) for m in mat[2].split("\t")[1:]]
        r = [r1, r2]
        M = {}
        for i, j in itt.product([0, 1], [0, 1]):
            lab = f"{c2[i]}{c1[j]}"
            M[lab] = r[i][j]
        return M

    # map
    mp = find_subsequent("Character to attribute mapping")[1:]
    mp = {m.split(":")[0]: m.split(":")[1].strip() for m in mp}
    info["map"] = mp

    # substitution rate
    sr = find_subsequent("Substitution rate")[0].split(" ")[-1]
    info["substitution_rate"] = float(sr)

    # equilibrium freq
    ef = find_subsequent("Equilibrium frequencies")[1:]
    ef = {mp[e.split(":")[0]]: float(e.split(":")[1].strip()) for e in ef}
    info["eq_freq"] = ef

    # symmetrized rates
    sr = find_subsequent("Symmetrized rates from")[1:]
    info["symm_rates"] = parse_matrix(mp, sr)

    # actual rates
    ar = find_subsequent("Actual rates from")[1:]
    info["act_rates"] = parse_matrix(mp, ar)

    return info


def parse_nexus(fname):
    pa_pattern = {}
    node_name = {}

    # load nexus tree
    tree = Phylo.read(fname, "nexus")

    # extract presence/absence pattern
    def parse_prop(n):
        prop = n.comment
        return prop.split('"')[-2]

    for n in tree.get_terminals():
        node_name[n] = n.name
        pa_pattern[n.name] = parse_prop(n)

    for n in tree.get_nonterminals():
        node_name[n] = n.confidence
        pa_pattern[n.confidence] = parse_prop(n)

    # find gain-loss events
    events = []

    def find_events(parent):
        if parent.is_terminal():
            return

        pn = node_name[parent]
        ps = pa_pattern[pn]
        for child in parent.clades:
            cn = node_name[child]
            cs = pa_pattern[cn]
            if cs != ps:
                if ps == "-":
                    events.append((cn, "gain"))
                elif cs == "-":
                    events.append((cn, "loss"))
                else:
                    events.append((cn, f"{cs}|{ps}"))
            find_events(child)

    find_events(tree.root)

    info = {"pa": pa_pattern, "events": events}
    return info


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


class TestClasses(unittest.TestCase):
    def test_node(self):
        n1 = Node("a", True)
        n2 = Node("a", False)

        self.assertEqual(n1, n1)
        self.assertNotEqual(n1, n2)
        self.assertEqual(hash(n1), hash(n1))
        self.assertNotEqual(hash(n1), hash(n2))
        self.assertEqual(n1.invert(), n2)
        self.assertEqual(n2.invert(), n1)
        self.assertEqual(hash(n1.invert()), hash(n2))
        self.assertNotEqual(n1, n1.invert())

    def test_path(self):
        n1 = Node("a", True)
        n2 = Node("b", True)
        n3 = Node("c", True)
        n4 = Node("d", True)

        p1 = Path([n1, n2, n3, n4])
        p2 = Path([n4.invert(), n3.invert(), n2.invert(), n1.invert()])

        self.assertEqual(p1, p1)
        self.assertEqual(p1, p2.invert())
        self.assertEqual(p2.invert(), p1)
        self.assertEqual(hash(p1), hash(p1))
        self.assertEqual(hash(p1), hash(p2.invert()))
        self.assertNotEqual(p1, p2)

        p3 = Path([n1, n2, n3])
        p3.add_right(n4)

        self.assertEqual(p1, p3)

        p4 = Path([n2, n3, n4])
        p4.add_left(n1)

        self.assertEqual(p1, p4)

    def test_edge(self):
        n1 = Node("a", True)
        n2 = Node("b", True)
        n3 = Node("c", True)

        e1 = Edge(n1, n2)
        e2 = Edge(n2.invert(), n1.invert())
        e3 = Edge(n1, n3)

        self.assertEqual(e1, e2)
        self.assertEqual(hash(e1), hash(e2))
        self.assertNotEqual(e1, e3)

    def test_junction(self):
        n1 = Node("a", True)
        n2 = Node("b", True)
        n3 = Node("c", True)
        n4 = Node("d", True)

        p1 = Path([n2, n3])
        p2 = Path([n3.invert(), n2.invert()])

        j1 = Junction(n1, p1, n4)
        j2 = Junction(n4.invert(), p2, n1.invert())

        self.assertEqual(j1, j1)
        self.assertEqual(j1, j2)
        self.assertEqual(hash(j1), hash(j1))
        self.assertEqual(hash(j1), hash(j2))

        self.assertTrue(j1.flanks_bid("a"))


if __name__ == "__main__":
    unittest.main()
