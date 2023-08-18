import argparse
import pathlib

import numpy as np
import pandas as pd

import pypangraph as pp
import utils as ut

from scipy.stats import entropy
from collections import defaultdict


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


def get_stats(pan_file):
    """
    Evaluate the following statistics:
    - n. blocks
    - n. nodes
    - n. path categories
    - majority category
    - whether there is only one path differing from the rest (sigleton)
    - how many paths only contain core blocks.
    - category entropy
    - min/max/mean path length
    - how many paths have duplicated blocks
    - left core length
    - right core length
    - average breakpoint entropy
    """
    # initialize stat dictionary
    stats = {"edge": pan_file.stem.split(".")[0]}

    # load pangraph
    pan = pp.Pangraph.load_json(pan_file)
    bdf = pan.to_blockstats_df()
    N = len(pan.strains())

    # number of blocks
    stats["n_blocks"] = len(bdf)
    stats["has_dupl"] = np.any(bdf["duplicated"].to_numpy())

    # number of categories
    paths = ut.pangraph_to_path_dict(pan)
    path_cat = path_categories(paths)
    stats["n_categories"] = len(path_cat)
    stats["majority_category"] = path_cat[0][0]
    stats["singleton"] = path_cat[0][0] == N - 1
    stats["cat_entropy"] = entropy([count for count, _, _ in path_cat])

    # path lengths
    Ls = []
    len_entr = 0
    n_all_cores = 0
    n_nodes = 0
    core_sides = {"left": [], "right": []}
    for count, nodes, isolates in path_cat:
        n_nodes += len(nodes)
        lengths = [bdf["len"][node.id] for node in nodes]
        Ls += [sum(lengths)] * count
        len_entr += entropy(lengths) * count
        core_sides["left"].append(nodes[0].id)
        core_sides["right"].append(nodes[-1].id)
        if np.all([bdf["core"][node.id] for node in nodes]):
            n_all_cores += count
    stats["n_nodes"] = n_nodes
    stats["min_length"] = min(Ls)
    stats["max_length"] = max(Ls)
    stats["mean_length"] = np.mean(Ls)
    stats["length_entropy"] = len_entr / N
    stats["n_all_cores"] = n_all_cores
    for side in ["left", "right"]:
        stats[f"core_{side}_length"] = None
        if not np.all(np.array(core_sides[side]) == core_sides[side][0]):
            continue
        if not bdf["core"][core_sides[side][0]]:
            continue
        stats[f"core_{side}_length"] = bdf["len"][core_sides[side][0]]

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
