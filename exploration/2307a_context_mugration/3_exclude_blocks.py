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


def describe_stratified_df(df, var):
    D = {}
    D["pangraph len"] = df.groupby(var).sum()["len"] / 1e6
    D["genome length"] = df.groupby(var).sum()["gen_len"] / 1e6 / len(strains)
    D["n. nodes"] = df.groupby(var).sum()["count"] / len(strains)
    D["n. blocks"] = df.groupby(var).count()["len"]

    for k, v in D.items():
        print(f"--- {k} ---")
        print(v, "\n--- percentage ---")
        print(v / v.sum(), "\n\n")


describe_stratified_df(bs_0, "duplicated")
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
describe_stratified_df(bs_1, "length category")

# %%
mask = bs_1["len"] >= ut.LEN_THR
bs_2 = bs_1[mask].copy().drop(columns=["length category"])
bs_2["core"].value_counts()

# %%
# remove synteny-break blocks

forbidden_blocks = np.loadtxt(ut.expl_fld / "forbidden_blocks.txt", dtype=str)
bs_2["synt. break"] = np.isin(bs_2.index, forbidden_blocks)

# %%
print("synteny breaks")
describe_stratified_df(bs_2, "synt. break")

print("synteny breaks accessory blocks")
describe_stratified_df(bs_2[~bs_2["core"]], "synt. break")

# %%

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

ax = axs[0]
sns.histplot(
    data=bs_2[~bs_2["core"]],
    x="n. strains",
    hue="synt. break",
    ax=ax,
    bins=np.arange(len(strains) + 1) + 0.5,
    element="step",
)
ax.set_xlabel("Number of strains")
ax.set_ylabel("Number of accessory  blocks")

ax = axs[1]
sns.histplot(
    data=bs_2[~bs_2["core"]],
    x="len",
    hue="synt. break",
    ax=ax,
    log_scale=True,
    element="step",
)
ax.set_xlabel("Block length (bp)")
ax.set_ylabel("Number of accessory blocks")

plt.tight_layout()
plt.savefig(ut.fig_fld / "3_synt_filter.png")
plt.show()

# %%
bs_3 = bs_2[~bs_2["synt. break"]].copy().drop(columns=["synt. break"])

retained_blocks = bs_3.index
np.savetxt(ut.expl_fld / "retained_blocks.txt", retained_blocks, fmt="%s")

# %%
