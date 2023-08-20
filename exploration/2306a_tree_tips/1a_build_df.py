# %%

import numpy as np
import pandas as pd
import pathlib
from Bio import Phylo
import pypangraph as pp

import utils as ut

from collections import defaultdict, Counter


svfld = pathlib.Path("data")
svfld.mkdir(exist_ok=True, parents=True)

# initialize dataframe
df = {}
# %%

prefix = "../../results/ST131/pangraph"
tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

# load tree and pangraph
tree = Phylo.read(tree_file, "newick")
pan = pp.Pangraph.load_json(pangraph_file)
# %%


# get length of terminal branches
def get_terminal_branch_lengths(tree):
    len_tips = {}
    for t in tree.get_terminals():
        len_tips[t.name] = t.branch_length
    return len_tips


df["branch_len"] = get_terminal_branch_lengths(tree)

# %%


def block_pa_measures(pan):
    """Extract block presence/absence measures for terminal branches from the pangraph"""

    bdf = pan.to_blockcount_df()  # dataframe with block-count per isolate
    is_dupl = pan.to_blockstats_df()["duplicated"].to_dict()
    len_dict = pan.to_blockstats_df()["len"].to_dict()

    # list of terminal blocks (present in only one isolate, possibly duplicated)
    mask = (bdf > 0).sum() == 1
    t_blocks = mask.index[mask].tolist()

    n_acc, n_dupl, L_acc, L_dupl = [defaultdict(int) for _ in range(4)]
    acc_list = defaultdict(list)  # list of acc blocks per isolate (non-duplicated)
    for i in pan.strains():
        for b in t_blocks:
            if bdf.loc[i, b] > 0:
                if is_dupl[b]:
                    n_dupl[i] += 1
                    L_dupl[i] += len_dict[b] * bdf.loc[i, b]
                else:
                    acc_list[i].append(b)
                    n_acc[i] += 1
                    L_acc[i] += len_dict[b]

    def n_block_groups(path, selected_blocks):
        """path is a list of blocks, with circular boundary conditions.
        This function counts group of adjacent selected blocks.
        """
        p = np.array(path)

        s = np.isin(p, selected_blocks).astype(int)
        breakpoints = np.abs(s - np.roll(s, 1))
        return np.sum(breakpoints) / 2

    n_events = {}
    for i, B in acc_list.items():
        p = pan.paths[i]
        # remove duplications from path
        bl_list = [b for b in p.block_ids if not is_dupl[b]]
        n_events[i] = n_block_groups(bl_list, B)

    return {
        "n_acc": n_acc,  # number of acc blocks
        "n_dupl": n_dupl,  # number of duplicated blocks
        "L_acc": L_acc,  # total length of acc blocks
        "L_dupl": L_dupl,  # total length of duplicated blocks
        "n_events": n_events,  # number of events
    }


df |= block_pa_measures(pan)

# %%


def path_to_dupl_junctions(path, is_dupl):
    """Parse a path and retrieve all duplicated junctions"""
    B = path.block_ids
    S = path.block_strands

    L = len(B)

    # make sure to start from non-duplicated block
    if is_dupl[B[0]]:
        # find last non-duplicated block:
        j = None
        for i in range(L - 1, 0, -1):
            if not is_dupl[B[i]]:
                j = i
                break
        # rotate path
        B = np.roll(B, -j)
        S = np.roll(S, -j)

    assert not is_dupl[B[0]], "First block is duplicated"

    J = []  # junctions
    nl, nr, p = None, None, None
    for i in range(L):
        b, s = B[i], S[i]

        if is_dupl[b]:
            if nl is None:
                # first duplicated block, initialize left node
                nl = ut.Node(B[i - 1], S[i - 1])
                p = ut.Path([ut.Node(b, s)])  # add first node to path
            else:
                p.add_right(ut.Node(b, s))  # continue path
        else:
            if nl is not None:
                # first non-duplicated block after duplication islands
                # terminate path and add junction
                nr = ut.Node(B[i], S[i])  # right node
                J.append(ut.Junction(nl, p, nr))  # create and add junction
                nl, nr, p = None, None, None  # reset current junction

    return J


def n_unique_duplicated_islands(pan):
    """Count the number of unique duplicated islands per isolate"""

    # list of all duplicated junctions per isolate
    is_dupl = pan.to_blockstats_df()["duplicated"].to_dict()
    junction_list = defaultdict(list)
    for i in pan.strains():
        junction_list[i] = path_to_dupl_junctions(pan.paths[i], is_dupl)

    # count all junctions
    j_counts = Counter([j for i in junction_list for j in junction_list[i]])

    # DEBUG
    # for i in j_counts:
    #     print(j_counts[i], i)

    # find all "terminal" unique junctions
    unique_j = [j for j in j_counts if j_counts[j] == 1]

    # find number of unique islands per isolate
    # using hashes makes comparison faster
    unique_hashes = [hash(j) for j in unique_j]
    unique_dupl_islands = {}
    for i in junction_list:
        hashes = [hash(j) for j in junction_list[i]]
        unique_dupl_islands[i] = np.isin(hashes, unique_hashes).sum()

    return unique_dupl_islands


df["n_terminal_dupl_islands"] = n_unique_duplicated_islands(pan)
# %%

# to pandas dataframe and save
df = pd.DataFrame(df).fillna(0)
df.to_csv(svfld / "tips_df.csv")
