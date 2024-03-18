import argparse
import pathlib
import subprocess
import tempfile

import pandas as pd

import utils as ut

from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--AB_states", type=str)
    parser.add_argument("--AB_nestedness", type=str)
    parser.add_argument("--joints_df", type=str)
    parser.add_argument("--tree", type=str)
    parser.add_argument("--out_nonsingleton_junct_df", type=str)
    parser.add_argument("--out_nonsingleton_branch_df", type=str)
    return parser.parse_args()


def mugration_inference(tree, AB_states):

    # create temporary folder
    with tempfile.TemporaryDirectory() as tmp_fld:

        # save tree copy
        tree_file = pathlib.Path(tmp_fld) / "core_tree.nwk"
        Phylo.write(
            tree,
            tree_file,
            "newick",
            format_branch_length="%1.5e",
        )

        # save AB states copy
        AB_states_file = pathlib.Path(tmp_fld) / "AB_states.csv"
        AB_states.to_csv(AB_states_file)

        tmp_path = pathlib.Path(tmp_fld)

        pa_inference = {}
        for i, state in enumerate(AB_states.columns):
            print(f"processing {state} \t - \t {i+1}/{len(AB_states.columns)}")
            out_dir = tmp_path / state

            cmd = f"""
            treetime mugration \
                --tree {tree_file} \
                --states {AB_states_file} \
                --attribute {state} \
                --outdir {out_dir} \
            """

            subprocess.run(cmd, shell=True)

            # output file names
            nexus_file = out_dir / "annotated_tree.nexus"
            rates_file = out_dir / "GTR.txt"

            # parse rates
            rates = ut.parse_gtr(rates_file)

            # parse inferred presence/absence pattern
            pa_pattern = ut.parse_nexus(nexus_file)

            # remove output directory
            subprocess.run(f"rm -r {out_dir}", shell=True)

            pa_inference[state] = {
                "pa_pattern": pa_pattern,
                "rates": rates,
            }

    return pa_inference


def initialize_branch_df(tree):
    branch_len = {b.name: b.branch_length for b in tree.find_clades()}
    branch_df = pd.DataFrame.from_dict(
        branch_len, orient="index", columns=["branch_length"]
    )
    branch_df["n_events"] = 0
    branch_df["terminal"] = False
    return branch_df


def analyze_mugration_output(tree, inf, AB_nestedness, AB_states):
    branch_df = initialize_branch_df(tree)
    event_info = {}
    for j, mug_inf in inf.items():
        info = {}
        nestedness = AB_nestedness["event_type"][j]
        events = mug_inf["pa_pattern"]["events"]
        info["n_events"] = len(events)
        info["n_minority"] = (AB_states[j] == "B").sum()
        info["gain"] = 0
        info["loss"] = 0
        info["undetermined"] = False
        for n, kind in events:
            branch_df.loc[n, "n_events"] += 1
            if nestedness == "A?B":
                info["undetermined"] = True
            elif nestedness == "A<B":
                if kind == "A|B":
                    info["loss"] += 1
                elif kind == "B|A":
                    info["gain"] += 1
            elif nestedness == "A>B":
                if kind == "A|B":
                    info["gain"] += 1
                elif kind == "B|A":
                    info["loss"] += 1

        event_info[j] = info
    idf = pd.DataFrame.from_dict(event_info, orient="index")
    return branch_df, idf


if __name__ == "__main__":

    args = parse_args()

    # select non-terminal events
    jdf = pd.read_csv(args.joints_df, index_col=0)
    mask = ~jdf["singleton"]
    jdf = jdf[mask]

    # load AB pattern
    AB_states = pd.read_csv(args.AB_states, index_col=0)
    AB_nestedness = pd.read_csv(args.AB_nestedness, index_col=0)

    # load core-genome tree
    tree = Phylo.read(args.tree, "newick")
    tree.ladderize()
    tree = ut.name_tree_nodes(tree)

    # perform mugration inference
    inf = mugration_inference(tree, AB_states)

    # analyze mugration output
    branch_df, idf = analyze_mugration_output(tree, inf, AB_nestedness, AB_states)

    # merge with other coldspot information
    jdf = pd.read_csv(args.joints_df, index_col=0)
    idf = jdf.merge(idf, left_index=True, right_index=True)

    # save dataframe
    idf.to_csv(args.out_nonsingleton_junct_df)
    branch_df.to_csv(args.out_nonsingleton_branch_df)
