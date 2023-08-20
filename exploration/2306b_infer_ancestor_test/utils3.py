import itertools as itt
import matplotlib.pyplot as plt

from Bio import Phylo


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
    event_type = {"P": "gain", "A": "loss"}

    def find_events(parent):
        if parent.is_terminal():
            return

        pn = node_name[parent]
        ps = pa_pattern[pn]
        for child in parent.clades:
            cn = node_name[child]
            cs = pa_pattern[cn]
            if cs != ps:
                events.append((cn, event_type[cs]))
            find_events(child)

    find_events(tree.root)

    info = {"pa": pa_pattern, "events": events}
    return info


def plot_tree_events(tree_file, pa_inference, bid):
    tree = Phylo.read(tree_file, "newick")

    info = pa_inference[bid]
    pa = info["pa_pattern"]["pa"]
    ev = info["pa_pattern"]["events"]

    def color_tree(node):
        # if node.is_terminal():
        #     if pa[node.name] == "P":
        #         node.color = "green"
        #     else:
        #         node.color = "red"
        # else:
        for nn, et in ev:
            if node.name == nn:
                node.color = "lime" if et == "gain" else "red"
        if node.color is None:
            node.color = "black"
        for c in node.clades:
            color_tree(c)

    def label_tree(node):
        if node.is_terminal():
            return node.name
        else:
            return ""

    def lab_colors(nn):
        if len(nn) == 0:
            return None
        if pa[nn] == "P":
            return "green"
        else:
            return "white"

    color_tree(tree.root)

    fig, ax = plt.subplots(1, 1, figsize=(10, 12))
    Phylo.draw(
        tree, label_func=label_tree, label_colors=lab_colors, axes=ax, do_show=False
    )
    plt.title(f"block - {bid}")
    return fig, ax
