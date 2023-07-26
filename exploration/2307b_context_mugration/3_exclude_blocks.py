# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import json

import utils as ut

# %%

# load pangraph
pan = ut.load_pangraph()
strains = pan.strains()

# block stats dataframe
bs_0 = pan.to_blockstats_df()

# %%

paths = ut.pangraph_to_path_dict(pan)

# %%

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

ax = axs[0]
sns.histplot(
    data=bs_0,
    x="n. strains",
    hue="duplicated",
    ax=ax,
    bins=np.arange(len(strains) + 1) + 0.5,
    element="step",
    # weights="len",
)
ax.set_xlabel("Number of strains")
ax.set_ylabel("Number of blocks")

ax = axs[1]
sns.histplot(
    data=bs_0,
    x="len",
    hue="duplicated",
    ax=ax,
    log_scale=True,
    element="step",
)
ax.set_xlabel("Block length (bp)")
ax.set_ylabel("Number of blocks")

plt.tight_layout()
plt.savefig(ut.fig_fld / "3_dupl_filter.png")
plt.show()

# %%

# how much duplicated sequence is removed ?
bs_0["gen_len"] = bs_0["len"] * bs_0["count"]
print("duplicated sequence removed:")
print(bs_0.groupby("duplicated").sum()["len"] / 1e6)
print(bs_0.groupby("duplicated").sum()["gen_len"] / 1e6 / len(strains))
print(bs_0.groupby("duplicated").sum()["count"])
print(bs_0.groupby("duplicated").count()["len"])


# %%
# remove sequence
mask = ~bs_0["duplicated"]
bs_1 = bs_0[mask].copy().drop(columns=["duplicated"])
bs_1["length category"] = (bs_1["len"] < ut.LEN_THR).replace(
    {True: "<500", False: ">=500"}
)
bs_1

# %%

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

ax = axs[0]
sns.histplot(
    data=bs_1,
    x="n. strains",
    hue="length category",
    ax=ax,
    bins=np.arange(len(strains) + 1) + 0.5,
    element="step",
)
ax.set_xlabel("Number of strains")
ax.set_ylabel("Number of blocks")

ax = axs[1]
sns.histplot(
    data=bs_1,
    x="len",
    hue="core",
    ax=ax,
    log_scale=True,
    element="step",
)
ax.axvline(ut.LEN_THR, color="k", linestyle="--")
ax.set_xlabel("Block length (bp)")
ax.set_ylabel("Number of blocks")

plt.tight_layout()
plt.savefig(ut.fig_fld / "3_len_filter.png")
plt.show()

# %%
print("short sequence removed:")
print(bs_1.groupby("length category").sum()["len"] / 1e6)
print(bs_1.groupby("length category").sum()["gen_len"] / 1e6 / len(strains))
print(bs_1.groupby("length category").sum()["count"])
print(bs_1.groupby("length category").count()["len"])

# %%


# def cleanup_paths(paths, keep_f):
#     """Given a function `keep_f` that takes block-ids as input,
#     remove nodes from paths that do not satisfy it."""
#     for iso, path in paths.items():
#         path = [node for node in path if keep_f(node.id)]
#         paths[iso] = path
#     return paths


# def merge_transitive_edges(paths, bdf):
#     pass

# %%
