from collections import defaultdict, Counter
import pathlib
import itertools as itt
from Bio import Phylo


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


class Junction:
    """A junction is a combination of a node (or path) flanked by two nodes, with reverse-complement simmetry"""

    def __init__(self, left: Node, center: Path, right: Node) -> None:
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

    def to_list(self):
        return [self.left.to_str_id(), self.center.to_list(), self.right.to_str_id()]

    def flanking_edge(self) -> Edge:
        return Edge(self.left, self.right)

    @staticmethod
    def from_list(t) -> "Junction":
        return Junction(
            Node.from_str_id(t[0]), Path.from_list(t[1]), Node.from_str_id(t[2])
        )


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


def path_junction_split(path, is_core):
    """Given a path and a boolean function of node ids, it splits the path in a set of
    Junctions, with flanking "core" blocks for which the condition is true."""
    junctions = []

    current = []
    left_node = None
    for node in path.nodes:
        if is_core(node.id):
            J = Junction(left_node, Path(current), node)
            junctions.append(J)
            left_node = node
            current = []
        else:
            current.append(node)

    # complete periodic boundary
    J = junctions[0]
    J.left = left_node
    J.center = Path(current + J.center.nodes)
    junctions[0] = J

    return junctions


def path_edge_count(paths):
    """Count internal edges of paths"""
    ct = Counter()
    for iso, p in paths.items():
        L = len(p.nodes)
        for i in range(L):
            e = Edge(p.nodes[i], p.nodes[(i + 1) % L])
            ct.update([e])
    return dict(ct)


def path_block_count(paths):
    """Count internal blocks of paths"""
    ct = Counter()
    for iso, p in paths.items():
        for node in p.nodes:
            ct.update([node.id])
    return dict(ct)


def find_mergers(paths):
    """Create a dictionary source -> sinks of block-ids to be merged"""
    edge_ct = path_edge_count(paths)
    block_ct = path_block_count(paths)

    mergers = {}
    for e, ec in edge_ct.items():
        bl, br = e.left.id, e.right.id
        if (ec == block_ct[bl]) and (ec == block_ct[br]):
            # merge
            if bl in mergers:
                if br in mergers:
                    source = mergers[br]
                    sink = mergers[bl]
                    for k in mergers:
                        if mergers[k] == source:
                            mergers[k] = sink
                else:
                    mergers[br] = mergers[bl]
            elif br in mergers:
                mergers[bl] = mergers[br]
            else:
                mergers[br] = bl
                mergers[bl] = bl
    return mergers


def create_pangraph_file_dict(pg_filenames):
    """
    Takes a list of pangraph junction filenames and returns a dictionary
    with the junction name as key and the pangraph filename as value.
    """
    pg_filenames = [pathlib.Path(p) for p in pg_filenames]
    return {p.stem: p for p in pg_filenames}


def name_tree_nodes(tree):
    """
    Assigns names to internal tree nodes. The names are assigned in a
    pre-order traversal of the tree.
    """
    for i, node in enumerate(tree.get_nonterminals(order="preorder")):
        node.name = f"int_node_{i}"
    return tree


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
                events.append((cn, f"{cs}|{ps}"))
            find_events(child)

    find_events(tree.root)

    info = {"pa": pa_pattern, "events": events}
    return info
