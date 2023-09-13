# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import numpy as np
import os
import json

import pypangraph as pp
import utils as ut
import mugration_utils as mu


from collections import defaultdict
from Bio import Phylo


fig_fld = pathlib.Path("figs")


fld = pathlib.Path("../../results/ST131")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

# %%
df = pd.read_csv(df_file, index_col=0)
Js = df.index.to_list()
# %%

tree = Phylo.read(tree_file, "newick")
tree.ladderize()
branch_len = {b.name: b.branch_length for b in tree.get_terminals()}
branch_len |= {b.name: b.branch_length for b in tree.get_nonterminals()}

# %%


def event_category(p1, p2):
    p1, p2 = set(p1), set(p2)
    if p1.issubset(p2):
        return "gain", list(p2 - p1)
    elif p2.issubset(p1):
        return "loss", list(p1 - p2)
    else:
        return "other", list((p1 | p2) - (p1 & p2))


def are_clade(tree, isolates):
    "check if the isolates form a clade in the tree"
    # get common ancestor
    leaves = {leaf.name: leaf for leaf in tree.get_terminals()}
    common_ancestor = tree.common_ancestor(*[leaves[i] for i in isolates])
    # check if terminal nodes are equal to clade
    clade = set([leaf.name for leaf in common_ancestor.get_terminals()])
    if clade == set(isolates):
        return common_ancestor.name
    else:
        return None


def find_coordinates(pan, block):
    for path in pan.paths:
        B = list(path.block_ids)
        if not block in B:
            continue
        idx = B.index(block)
        strand = path.block_strands[idx]
        beg, end = path.block_positions[idx : idx + 2]
        iso = path.name
        return beg, end, iso, strand


events_df = {}
for j in Js:
    info = {}
    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    paths = ut.pangraph_to_path_dict(pan)
    path_cat = ut.path_categories(paths)

    if not len(path_cat) == 2:
        continue

    n1, p1, i1 = path_cat[0]
    n2, p2, i2 = path_cat[1]

    iso = None
    if n2 == 1:
        info["type"] = "terminal"
        iso = i2[0]
    else:
        info["type"] = "internal"
        iso = are_clade(tree, i2)

    if not (iso is None):
        info["isolate"] = iso
        info["branch len"] = branch_len[iso]

    cat, blocks = event_category(p1, p2)
    info["event category"] = cat
    if len(blocks) == 1:
        bid = blocks[0].id
        info["block"] = bid
        info["seq"] = pan.blocks[bid].sequence
        bdf = pan.to_blockstats_df()
        info["block len"] = bdf.loc[bid, "len"]

        beg, end, loc_iso, strand = find_coordinates(pan, bid)
        info["loc_beg"] = beg
        info["loc_end"] = end
        info["loc_iso"] = loc_iso
        info["loc_strand"] = strand

    events_df[j] = info

events_df = pd.DataFrame.from_dict(events_df, orient="index")
events_df.to_csv("data/twocat_events.csv")
events_df
# %%

sns.histplot(events_df, x="event category", hue="type", multiple="stack")
plt.show()

# %%
sns.histplot(
    events_df,
    x="block len",
    hue="event category",
    log_scale=True,
    stat="probability",
    common_norm=False,
    common_bins=True,
    # multiple="stack",
)
plt.show()

# %%

# find annotations
from Bio import SeqIO
import json

gbk_fld = pathlib.Path("../../data/gbk")

offsets_file = "../../results/ST131/backbone_joints/asm20-100-5/joints_pos.json"
with open(offsets_file, "r") as f:
    offsets = json.load(f)


def parse_features(gbk_file, ws, we):
    assert ws < we, "start must be smaller than end"
    features_in_window = []

    for record in SeqIO.parse(gbk_file, "genbank"):
        L = len(record)
        if we > L:
            print("### wrapping ###")
            we = we % L
            if ws > L:
                ws = ws % L
                return parse_features(gbk_file, ws, we)
            else:
                return parse_features(gbk_file, ws, L) + parse_features(gbk_file, 0, we)

        for feature in record.features:
            if feature.type == "CDS":
                fs = feature.location.start.position
                fe = feature.location.end.position

                cs = max([fs, ws])
                ce = min([fe, we])

                if not (cs < ce):
                    continue

                left_within = ws < fs
                right_within = fe < we
                if left_within and right_within:
                    cat = "full"
                elif not (left_within or right_within):
                    cat = "within"
                else:
                    cat = "partial"

                overlap = ce - cs + 1
                # f_seq = str(feature.extract(record.seq))
                features_in_window.append((feature, overlap, cat))
        break
    # export sequence, feature size and window overlap
    return features_in_window


def parse_qualifiers(qual):
    dq = dict(qual)
    for k in [
        "translation",
        "transl_except",
        "note",
        "old_locus_tag",
        "inference",
        "ribosomal_slippage",
        "codon_start",
        "transl_table",
    ]:
        if k in dq:
            del dq[k]
    for k in dq:
        if isinstance(dq[k], list):
            if len(dq[k]) == 1:
                dq[k] = dq[k][0]
    return dq


ft_df = []
for j, row in events_df.iterrows():
    print(f"processing {j}")
    if not pd.isna(row["loc_beg"]):
        iso = row["loc_iso"]
        gbk_file = gbk_fld / f"{iso}.gbk"
        ws, we = row["loc_beg"], row["loc_end"]
        os, oe, o_strand = offsets[j][iso]
        ws += os
        we += os
        features = parse_features(gbk_file, ws, we)

        for feat, overlap, cat in features:
            entry = {
                "overlap": overlap,
                "cat": cat,
                "joint": j,
                "iso": iso,
                "start": int(feat.location.start),
                "end": int(feat.location.end),
                "feature_len": len(feat),
                "window_size": we - ws,
            }
            entry |= parse_qualifiers(feat.qualifiers)
            ft_df.append(entry)
ft_df = pd.DataFrame(ft_df)
ft_df.to_csv("data/annotations.csv")
ft_df

# %%
