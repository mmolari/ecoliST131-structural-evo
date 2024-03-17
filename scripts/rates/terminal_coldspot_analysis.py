import pandas as pd
import argparse
import pathlib

import pypangraph as pp
import utils as ut

from Bio import Phylo


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coldspots_df", type=str)
    parser.add_argument("--tree", type=str)
    parser.add_argument("--junction_pangraphs", type=str, nargs="+")
    parser.add_argument("--out_coldspot_df", type=str)
    parser.add_argument("--out_events_df", type=str)

    return parser.parse_args()


def count_events(path_cats):
    """
    Takes as input paths that:
    - have two categories
    - the minority category is present in only one isolate
    Decides whether the minority category was subject to a gain or a loss.
    """
    assert len(path_cats) == 2
    p1, p2 = path_cats
    c1, n1, i1 = p1
    c2, n2, i2 = p2
    assert c2 == 1, f"{c1}, {c2}"
    assert n1[0] == n2[0], f"{n1}, {n2}"
    assert n1[-1] == n2[-1], f"{n1}, {n2}"
    assert len(i2) == 1
    i2 = i2[0]
    n1, n2 = set(n1), set(n2)

    if n1.issubset(n2):
        return i2, "gain", n2 - n1
    elif n2.issubset(n1):
        return i2, "loss", n1 - n2
    else:
        return i2, "other", (n1 | n2) - (n1 & n2)


def create_pangraph_file_dict(pg_filenames):
    """
    Takes a list of pangraph junction filenames and returns a dictionary
    with the junction name as key and the pangraph filename as value.
    """
    pg_filenames = [pathlib.Path(p) for p in pg_filenames]
    return {p.stem: p for p in pg_filenames}


if __name__ == "__main__":
    args = parse_args()

    # select only singleton coldspots
    cdf = pd.read_csv(args.coldspots_df, index_col=0)
    terminal_mask = cdf["singleton"]
    cdf = cdf[terminal_mask]

    # parse tree
    tree = Phylo.read(args.tree, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    terminal_len = {b.name: b.branch_length for b in tree.get_terminals()}

    # pangraph filenames
    pangraph_files = create_pangraph_file_dict(args.junction_pangraphs)

    edf = pd.DataFrame.from_dict(
        terminal_len, orient="index", columns=["branch_length"]
    )
    edf["gain"] = 0
    edf["loss"] = 0
    edf["other"] = 0

    # add columns to coldspots df
    cdf["event_type"] = None
    cdf["event_iso"] = None

    # for every coldspot
    for j in cdf.index.to_list():

        # load pangraph
        pan = pp.Pangraph.load_json(pangraph_files[j])

        # extract paths and categories
        paths = ut.pangraph_to_path_dict(pan)
        path_cat = ut.path_categories(paths)

        # categorize as gain/loss event and isolate
        iso, tp, bs = count_events(path_cat)

        # append results
        edf.loc[iso, tp] += 1
        cdf.loc[j, "event_type"] = tp
        cdf.loc[j, "event_iso"] = iso

    edf.to_csv(args.out_coldspot_df)
    cdf.to_csv(args.out_events_df)
