# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import pathlib
import json
from Bio import Phylo, SeqIO

# fig_fld = pathlib.Path("figs/f04")
# fig_fld.mkdir(exist_ok=True, parents=True)

data_fld = pathlib.Path("data/f04")
data_fld.mkdir(exist_ok=True, parents=True)

dset = "ST131_ABC"
fld = pathlib.Path(f"../../results/{dset}")

# load annotations
fld_data = fld / "annotations/junct_pos_asm20-100-5"
is_df = pd.read_csv(fld_data / "ISEScan_real.csv", index_col=0)

# load info
is_info = pd.read_csv(fld / "annotations/isescan/is_summary.tsv", sep="\t")
is_info["id"] = (
    is_info["seqID"] + "|" + is_info["family"] + "|" + is_info.index.astype(str)
)
is_info = is_info.set_index("id", verify_integrity=True)


# load joint coordinates dictionary
pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
with open(pos_file) as f:
    jp = json.load(f)

# load junction info dataframes
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_j = pd.read_csv(df_file, index_col=0)
df_j["delta_len"] = df_j["max_length"] - df_j["min_length"]

# load tree
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]

# genome lengths
iso_L_file = fld / "pangraph/genome_lengths.csv"
iso_L = pd.read_csv(iso_L_file, index_col=0)["length"].to_dict()

# %%

# add number of IS
df_j["n_IS"] = is_df["junction"].value_counts()
df_j["n_IS"] = df_j["n_IS"].fillna(0)
df_j["n_IS"] = df_j["n_IS"].astype(int)

mask = df_j["n_IS"] > 0
mask &= df_j["n_categories"] == 2
mask &= df_j["n_iso"] == 222
sdf = df_j[mask]
Js = sdf.index.to_list()
sdf
# 280 backbone junctions with 2 categories and at least 1 IS (there were 3 extra common)


# %%

# n_IS  majority_category  count  compatible with gain/loss
# 1     221                  217  x
# 2     220                   21  x
#       221                   13
# 3     221                    6
#       219                    3  x
# 4     218                    3  x
# 21    201                    2  x
# 6     216                    2  x
# 8     214                    2  x
# 168   168                    1  x (loss)
# 162   162                    1  x (loss)
# 52    170                    1  x
# 23    199                    1  x
# 7     221                    1
# 17    221                    1
# 9     221                    1
# 5     221                    1
# 4     220                    1
# 3     220                    1
# 223   221                    1

# (254 seem compatible)

# for each junction:
# - consistency check
#   - check that the isolates where the accessory region is empty are the ones without IS
#   - and the ones with accessory region are the ones with IS
# - extract all annotations in empty region
#   - check annotation consistency
# - save:
#   - junction id
#   - number of IS/empty (fitness?)
#   - whether a gene was interrupted
#   - if so, which gene

# %%


def empty_full(pos, iso_L):
    empty, full, empty_pos = [], [], {}
    for iso in pos:
        cb, ab, ae, ce, strand = pos[iso]
        if (ae - ab) % iso_L[iso] == 1:
            empty.append(iso)
            empty_pos[iso] = ab
        else:
            print("nonempy ->", iso, " | accessory L =", (ae - ab) % iso_L[iso])
            full.append(iso)
    return empty, full, empty_pos


def has_IS(j, is_df):
    mask = is_df["junction"] == j
    return is_df[mask]["iso"].unique()


def extract_gbk_annotation(gbk_file, pos):
    record = SeqIO.read(gbk_file, "genbank")
    hits = []
    for feature in record.features:
        if feature.type == "source":
            continue
        s, e = int(feature.location.start), int(feature.location.end)
        if s > e:
            s, e = e, s
        L = np.abs(e - s)
        isin = False
        if L < 1e5:
            isin = (pos > s) and (pos < e)
        else:
            isin = pos in feature.location
        if isin:
            hits.append(
                {
                    "type": feature.type,
                    "start": s,
                    "end": e,
                    "strand": feature.location.strand,
                    "iso": record.id,
                    "pos": pos,
                    "gene": feature.qualifiers.get("gene", [None])[0],
                    "product": feature.qualifiers.get("product", [None])[0],
                    "GO_function": feature.qualifiers.get("GO_function", [None])[0],
                    "GO_process": feature.qualifiers.get("GO_process", [None])[0],
                }
            )
    return hits


import sys
import time


def process_junction(j, sdf=sdf, jp=jp, iso_L=iso_L, is_df=is_df):

    t0 = time.time()

    n_IS = sdf.loc[j, "n_IS"]
    maj_cat = sdf.loc[j, "majority_category"]
    print(f"{j=}, {n_IS=}, {maj_cat=}")
    pos = jp[j]

    # check empty junction
    empty, full, empty_pos = empty_full(pos, iso_L)

    # check has IS
    is_isolates = has_IS(j, is_df)

    s_res = {
        "junction": j,
        "n_empty": len(empty),
        "n_full": len(full),
        "n_IS": n_IS,
        "n_IS_isolates": len(is_isolates),
    }

    # check consistency
    if set(full) == set(is_isolates):
        s_res["nonempty_is_IS"] = True
        print("nonempty == IS clones --- check passed")
    else:
        s_res["nonempty_is_IS"] = False
        print("nonempty != IS clones --- check failed /// skipping\n")
        return s_res

    # check gene interruption)
    hits = []
    for iso, p in empty_pos.items():
        gbk_file = f"../../data/gbk/{iso}.gbk"
        hits += extract_gbk_annotation(gbk_file, p)
    hits = pd.DataFrame(hits)

    # if no hits
    if len(hits) == 0:
        print("no hits found /// skipping\n")
        return s_res
    print(f"hits found - {len(hits)}")

    for tp, n in hits["type"].value_counts().items():
        s_res[f"n_hits_{tp}"] = n
        print(f"n. hits {tp}: {n}")
    hits.to_csv(data_fld / f"hits_{j}.csv")

    cds_hits = hits[hits["type"] == "CDS"]
    # check annotation consistency
    if cds_hits["product"].nunique() == 1:
        tp = "consistent"
        print("annotation consistency check passed")
    else:
        tp = "inconsistent"
        print("annotation consistency check failed")
    s_res["n_products"] = cds_hits["product"].nunique()
    s_res["annotation_consistency"] = tp
    s_res["annotations"] = "|".join(cds_hits["product"].unique())

    print(f"annotations: {cds_hits['product'].unique()}\n")
    tf = time.time()
    print(f"Elapsed time: {tf - t0:.2f} s\n")
    return s_res


# def progressbar(it, prefix="", size=60, out=sys.stdout):  # Python3.6+
#     count = len(it)
#     start = time.time()

#     def show(j):
#         x = int(size * j / count)
#         remaining = ((time.time() - start) / j) * (count - j)

#         mins, sec = divmod(remaining, 60)
#         time_str = f"{int(mins):02}:{sec:05.2f}"

#         print(
#             f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count} Est wait {time_str}",
#             end="\r",
#             file=out,
#             flush=True,
#         )

#     for i, item in enumerate(it):
#         yield item
#         show(i + 1)
#     print("\n", flush=True, file=out)


# summary_df = []
# for j in progressbar(Js):
#     res = process_junction(j, sdf, jp, iso_L, is_df)
#     summary_df.append(res)

import multiprocessing

with multiprocessing.Pool(20) as pool:
    summary_df = pool.map(process_junction, Js)

summary_df = pd.DataFrame(summary_df)
summary_df.to_csv(data_fld / "summary.csv")

# %%
