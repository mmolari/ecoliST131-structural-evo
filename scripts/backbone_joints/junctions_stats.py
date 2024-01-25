import argparse
import pathlib

import numpy as np
import pandas as pd

import pypangraph as pp
import utils as ut

from scipy.stats import entropy


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Given a set of junction pangraphs, compute statistics for the junctions and
        saves them in a dataframe.
        """
    )
    parser.add_argument("--junct_pangraphs", type=str, nargs="+")
    parser.add_argument("--df_csv", type=str)
    return parser.parse_args()


def get_stats(pan_file):
    """
    Evaluate the following statistics:
    - n. blocks
    - n. nodes
    - n. path categories
    - n. isolates
    - majority category
    - whether there is only one path differing from the rest (sigleton)
    - how many paths only contain core blocks.
    - category entropy
    - min/max/mean path length
    - how many paths have duplicated blocks
    - whether the joint is transitive (only one path category)
    - average accessory non-empty length
    - left core length
    - right core length
    """
    # initialize stat dictionary
    stats = {"edge": pan_file.stem.split(".")[0]}

    # load pangraph
    pan = pp.Pangraph.load_json(pan_file)
    bdf = pan.to_blockstats_df()
    N = len(pan.strains())
    stats["n_iso"] = N

    # number of blocks
    stats["n_blocks"] = len(bdf)
    stats["has_dupl"] = np.any(bdf["duplicated"].to_numpy())

    # number of categories
    paths = ut.pangraph_to_path_dict(pan)
    path_cat = ut.path_categories(paths)
    stats["n_categories"] = len(path_cat)
    stats["majority_category"] = path_cat[0][0]
    # singleton: only one path differing from the rest
    stats["singleton"] = path_cat[0][0] == N - 1
    stats["cat_entropy"] = entropy([count for count, _, _ in path_cat])

    # path lengths
    Ls = []
    len_entr = 0
    n_all_cores = 0
    n_nodes = 0
    core_sides = {"left": [], "right": []}
    # average length of accessory blocks
    for count, nodes, isolates in path_cat:
        n_nodes += len(nodes)
        lengths = [bdf["len"][node.id] for node in nodes]
        Ls += [sum(lengths)] * count
        core_sides["left"].append(nodes[0].id)
        core_sides["right"].append(nodes[-1].id)
        if np.all([bdf["core"][node.id] for node in nodes]):
            n_all_cores += count
    stats["n_nodes"] = n_nodes
    stats["min_length"] = min(Ls)
    stats["max_length"] = max(Ls)
    stats["mean_length"] = np.mean(Ls)
    stats["n_all_cores"] = n_all_cores
    for side in ["left", "right"]:
        stats[f"core_{side}_length"] = None
        if not np.all(np.array(core_sides[side]) == core_sides[side][0]):
            continue
        if not bdf["core"][core_sides[side][0]]:
            continue
        stats[f"core_{side}_length"] = bdf["len"][core_sides[side][0]]

    # whether the junction only contains a single core block, no variation
    stats["transitive"] = stats["n_categories"] == 1

    # average length of accessory blocks when non-empty
    # and average frequency of only core blocks.
    acc_len = 0
    n_with_acc = 0
    for count, nodes, isolates in path_cat:
        acc_lengths = [
            bdf["len"][node.id] for node in nodes if not bdf["core"][node.id]
        ]
        if len(acc_lengths) > 0:
            acc_len += sum(acc_lengths) * count
            n_with_acc += count
    if n_with_acc > 0:
        stats["nonempty_acc_len"] = acc_len / n_with_acc
    else:
        stats["nonempty_acc_len"] = None
    stats["nonempty_freq"] = n_with_acc / N

    return stats


if __name__ == "__main__":
    args = parse_args()

    # capture files
    pan_files = [pathlib.Path(f) for f in args.junct_pangraphs]

    # extract stats
    df = []
    for pan_file in pan_files:
        stats = get_stats(pan_file)
        df.append(stats)
    df = pd.DataFrame(df)

    # save to csv
    df.to_csv(args.df_csv, index=False)
