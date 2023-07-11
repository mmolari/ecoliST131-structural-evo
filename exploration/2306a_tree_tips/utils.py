import unittest


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


class Junction:
    """A junction is a combination of a path flanked by two nodes, with reverse-complement simmetry"""

    def __init__(self, left: Node, center: Path, right: Node) -> None:
        self.left = left
        self.center = center
        self.right = right

    def invert(self) -> "Junction":
        return Junction(self.right.invert(), self.center.invert(), self.left.invert())

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


if __name__ == "__main__":
    unittest.main()