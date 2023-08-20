# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %%
df = pd.read_csv(
    "../../results/ST131/backbone_joints/asm20-100-5/junctions_stats.csv",
    index_col=0,
)

N = 84

df["delta_len"] = df["max_length"] - df["min_length"]
df["delta_len_cat"] = pd.cut(df["delta_len"], [0, 100, 1000, 10000, 100000])

df["n_blocks_cat"] = pd.cut(df["n_blocks"], [0, 10, 100])

df["type"] = "complex"
df.loc[df["n_blocks"] == 1, "type"] = "no variation"
df.loc[df["singleton"], "type"] = "singleton (exchange)"
df.loc[df["singleton"] & (df["n_all_cores"] == (N - 1)), "type"] = "singleton (gain)"
df.loc[df["singleton"] & (df["n_all_cores"] == 1), "type"] = "singleton (loss)"

df["type"] = pd.Categorical(
    df["type"],
    categories=[
        "no variation",
        "singleton (exchange)",
        "singleton (gain)",
        "singleton (loss)",
        "complex",
    ],
)


sns.histplot(
    data=df,
    x="delta_len",
    hue="type",
    element="step",
    # multiple="stack",
    bins=[-0.5, 0.5, 1.5] + list(np.logspace(1, 5, 20)),
    common_norm=False,
    stat="probability",
)
plt.xscale("symlog", linthresh=10)
plt.xlim(left=-0.5)
plt.show()

# columns = [
#     "n_blocks",
#     "has_dupl",
#     "n_categories",
#     "majority_category",
#     "singleton",
#     "cat_entropy",
#     "n_nodes",
#     "min_length",
#     "max_length",
#     "mean_length",
#     "length_entropy",
#     "n_all_cores",
#     "core_left_length",
#     "core_right_length",
#     "delta_len",
# ]
# %%
df.value_counts(["type"])

# %%
sdf = df[df["singleton"]]
sns.histplot(data=sdf, x="n_all_cores", element="step")
plt.yscale("log")
plt.show()

sdf["n_all_cores"].value_counts()

#  'HBNSCZZTNM_f__JEWZJAMVQY_f' -> bigger (gain)
#  'GSDMALBSCR_f__NQAQWJWJBU_f' -> same
#  'NGGVZDRMNK_r__RVFLSKSCEH_f' -> smaller (loss)
#  'KIHCQKPRTH_r__QIAUDBKRYS_r' -> smaller (loss)
#  'FOYTPQJHOG_r__OXWYRVTJCX_f' -> smaller (loss, duplicated)
#  'MZQNCXKOQO_r__PNEDIIFRCD_r' -> weird, moved to the left.
# %%
sdf = df[(~df["singleton"]) & (df["n_categories"] > 1)]

sns.histplot(data=sdf, x="n_categories", element="step")
plt.yscale("log")
plt.show()

sns.histplot(data=sdf, x="majority_category", element="step")
plt.yscale("log")
plt.show()

# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
sns.scatterplot(
    data=sdf,
    x="n_categories",
    y="majority_category",
    hue="delta_len_cat",
    # hue="has_dupl",
    alpha=0.3,
    ax=ax,
)
plt.legend(loc="upper right", title="(max - min) lenght")
plt.xlabel("Number of path categories")
plt.ylabel("Majority path category")
plt.grid(True, alpha=0.5, ls=":")
ax.set_xlim(left=0)
sns.despine()
plt.tight_layout()
plt.savefig("figs/scr_5/delta_len.png", facecolor="white")
plt.show()

# %%
df[df["majority_category"] < 60]
# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
sns.scatterplot(
    data=sdf,
    x="n_categories",
    y="majority_category",
    hue="has_dupl",
    # hue="has_dupl",
    alpha=0.3,
    ax=ax,
)
plt.legend(loc="upper right", title="has duplicated blocks")
plt.xlabel("Number of path categories")
plt.ylabel("Majority path category")
plt.grid(True, alpha=0.5, ls=":")
ax.set_xlim(left=0)
sns.despine()
plt.tight_layout()
plt.savefig("figs/scr_5/has_dupl.png", facecolor="white")
plt.show()

# %%

sns.scatterplot(
    data=df,
    x="n_categories",
    y="n_blocks",
    hue="delta_len_cat",
    alpha=0.5,
)
plt.legend(loc="lower right", title="(max - min) lenght")
plt.xlabel("Number of path categories")
plt.ylabel("Number of blocks")
plt.yscale("log")
plt.grid(True, alpha=0.5, ls=":")
sns.despine()
plt.tight_layout()
plt.savefig("figs/scr_5/n_blocks.png", facecolor="white")
plt.show()

# %%
