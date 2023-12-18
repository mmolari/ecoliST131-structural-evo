import os
import json
import pathlib
import itertools as itt
import pandas as pd
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


def mugration_inference(tree_file, states_file, svfld, inference_svfile):
    pa_inference = {}

    svfld = pathlib.Path(svfld)

    state_list = pd.read_csv(states_file, index_col=0).columns.to_list()

    for i, state in enumerate(state_list):
        print(f"processing {state} \t - \t {i+1}/{len(state_list)}")
        out_dir = svfld / state

        cmd = f"""
        treetime mugration \
            --tree {tree_file} \
            --states {states_file} \
            --attribute {state} \
            --outdir {out_dir} \
        """

        os.system(cmd)

        # output file names
        nexus_file = out_dir / "annotated_tree.nexus"
        rates_file = out_dir / "GTR.txt"

        # parse rates
        rates = parse_gtr(rates_file)

        # parse inferred presence/absence pattern
        pa_pattern = parse_nexus(nexus_file)

        os.system(f"rm -r {out_dir}")

        pa_inference[state] = {
            "pa_pattern": pa_pattern,
            "rates": rates,
        }

    # write to json
    with open(inference_svfile, "w") as f:
        json.dump(pa_inference, f)
